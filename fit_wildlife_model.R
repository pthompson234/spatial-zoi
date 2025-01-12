source("FUNCTIONS.R")

require(TMB)
require(terra)
require(data.table)
require(dplyr)
require(readr)
require(magrittr)
require(minqa)
require(optimx)
require(tmbstan)

SCALE_ELEV = 1000 # m to km
SCALE_DTWN = 1000 # m to km
SCALE_ALPHA_DIST = 1000 # m to km
MAX_DIST_TO_KEEP = 12000 # in meters; how far away do we stop including disturbance values? A value of 12000 is very conservative and takes a while, but may be necessary for animals with low tolerance for humans
D_ROUND = 30 # in meters; how do we group together the distances?

SEASON_PARS = c(106.28977, 307.90013, 18.58934, 15.14239) # Seasonal activity parameters (see Supplementary Material for more information)

detections = read_csv("grizzly_detections_SAMPLE.csv")
sites = read_csv("grizzly_camera_effort_SAMPLE.csv")
all_site_data = read_csv("site_data_SAMPLE.csv")

## FILTER OUT SITES THAT ARE NOT IN STUDY AREA
all_site_data = all_site_data %>% dplyr::filter(in_study == 1, uptime_days > 1)
detections = detections %>% dplyr::filter(site_name %in% all_site_data$site_name)
sites = sites %>% dplyr::filter(site_name %in% all_site_data$site_name)

# Add covariates for everything (scaling / squaring / etc)
detections = detections %>% mutate(intercept = 1,
                                   elevation_km = elevation / SCALE_ELEV, 
                                   slope2 = slope^2, 
                                   cos_aspect = cos(aspect_radians), 
                                   sin_aspect = sin(aspect_radians), 
                                   dtown_km = dtown_m / SCALE_DTWN, 
                                   droad_km = droad_m / SCALE_DTWN, 
                                   snow_annual_2 = snow_annual ^ 2) %>% 
  mutate(elevation2_km = elevation_km ^ 2)

sites = sites %>% mutate(intercept = 1,
                         elevation_km = elevation / SCALE_ELEV, 
                         slope2 = slope^2, 
                         cos_aspect = cos(aspect_radians), 
                         sin_aspect = sin(aspect_radians), 
                         dtown_km = dtown_m / SCALE_DTWN,
                         droad_km = droad_m / SCALE_DTWN, 
                         snow_annual_2 = snow_annual ^ 2) %>%
  mutate(elevation2_km = elevation_km ^ 2)

site_histories = sites[, c("DateStart", "DateLastWorking")] %>% as.matrix

all_site_data = all_site_data %>% mutate(intercept = 1,
                                         elevation_km = elevation / SCALE_ELEV, 
                                         slope2 = slope^2, 
                                         cos_aspect = cos(aspect_radians), 
                                         sin_aspect = sin(aspect_radians), 
                                         dtown_km = dtown_m / SCALE_DTWN,
                                         droad_km = droad_m / SCALE_DTWN, 
                                         snow_annual_2 = snow_annual ^ 2) %>%
  mutate(elevation2_km = elevation_km ^ 2)

t_weights_log_sum = sum(log(seasonality_func(detections$t, SEASON_PARS, mod = TRUE)))
site_weights = int_seasonality_func(sites$DateStart, sites$DateLastWorking, SEASON_PARS, mod = TRUE)

covariates_E = c("intercept", "elevation_km", "elevation2_km", "slope", "slope2", "cos_aspect", "sin_aspect", "snow_annual", "snow_annual_2", "barren", "herbaceous", "conif_open", "shrub")

t_covariates_sum = colSums(detections[, covariates_E])
site_covariates = as.matrix(sites[, covariates_E])

# get distances between sites by converting to SpatVector
all_site_data_sp = vect(all_site_data, geom = c("easting", "northing"), crs = "+init=epsg:26911")

# This is the actual output generated by the full dataset. For use in the wildlife models
human_disturbance = rast("humantraffic.tif")
names(human_disturbance) = "humans_per_day"

message("Getting zonal statistics for human disturbance...")

# This takes a bit of time, which can be tracked using the progress bar.
pb = txtProgressBar(style = 3)
all_site_zonal_values = lapply(1:nrow(all_site_data_sp), function(row) {

  this_site = all_site_data_sp[row, ]

  cropped_raster = crop(human_disturbance, ext(this_site) + MAX_DIST_TO_KEEP)

  d_this_site = distance(cropped_raster, this_site)
  d_this_site[d_this_site > MAX_DIST_TO_KEEP] = NA

  zobj_orig = zonal(cropped_raster, d_this_site, fun = "sum") %>% as.data.frame
  setDT(zobj_orig)
  zobj_orig = zobj_orig[, dist_floor := ceiling(humans_per_day / D_ROUND) * D_ROUND]
  zobj_orig$dist_floor[1] = 0 # the first value always should represent the raster cell we are in
  zobj = zobj_orig[, .(sum(humans_per_day.1)), by = .(dist_floor)]
  zobj_dists = zobj_orig[, .(mean(humans_per_day)), by = .(dist_floor)]
  zobj_dists[1, ] = data.table(0, all_site_data$humans_per_day[row])

  setTxtProgressBar(pb, value = row / nrow(all_site_data_sp))

  cbind(zobj_dists[, 2], zobj[, 2])
})
close(pb)

site_zonal_values_matrix = do.call("cbind", all_site_zonal_values) %>% as.matrix

all_site_distances = site_zonal_values_matrix[, seq(1, ncol(site_zonal_values_matrix), by = 2)]
all_site_covariates = site_zonal_values_matrix[, seq(2, ncol(site_zonal_values_matrix), by = 2)]

# get the row index of our all_site_data data.frame for each detection and site history
t_site_indices = sapply(1:nrow(detections), function(row) which(all_site_data$site_name == detections$site_name[row]))
site_site_indices = sapply(1:nrow(sites), function(row) which(all_site_data$site_name == sites$site_name[row]))

t_distances = all_site_distances[, t_site_indices] %>% t
site_distances = all_site_distances[, site_site_indices] %>% t

# turn meters into kilometers (or whatever unit is desired)
t_distances = t_distances / SCALE_ALPHA_DIST
site_distances = site_distances / SCALE_ALPHA_DIST

t_spcov_values = all_site_covariates[, t_site_indices] %>% t
site_spcov_values = all_site_covariates[, site_site_indices] %>% t

# for debugging purposes - a test
SPCOV_SCALE = 100

t_spcov_values = t_spcov_values / SPCOV_SCALE
site_spcov_values = site_spcov_values / SPCOV_SCALE

# Correct zero's (or negative numbers) because they are not suitable for this model! There are basically none of these in our case, here, and because there is (in reality) tiny amounts of human disturbance just about everywhere in AB, I don't think this is very consequential.
MIN_SPCOV_VALUE = 0.0001

t_spcov_values[t_spcov_values <= 0] = MIN_SPCOV_VALUE
site_spcov_values[site_spcov_values <= 0] = MIN_SPCOV_VALUE

# convert dimensions so it's a 3-D array even if there is only one covariate
if (length(dim(t_spcov_values)) == 2) dim(t_spcov_values) = c(dim(t_spcov_values), 1)
if (length(dim(site_spcov_values)) == 2) dim(site_spcov_values) = c(dim(site_spcov_values), 1)

covariates_A = c("intercept") # overwrite what used to be the intercept

t_alcov_values = as.matrix(detections[, covariates_A, drop = FALSE])
site_alcov_values = as.matrix(sites[, covariates_A, drop = FALSE])

model_fit_lamE_betH_alpA = fit_ID_TMB(list(t_weights_log_sum = t_weights_log_sum,
                                           t_covariates_sum = t_covariates_sum,
                                           t_distances = t_distances,
                                           t_spcov_values = t_spcov_values,
                                           t_alcov_values = t_alcov_values,
                                           site_weights = site_weights,
                                           site_covariates = site_covariates,
                                           site_distances = site_distances,
                                           site_spcov_values = site_spcov_values,
                                           site_alcov_values = site_alcov_values),
                                      method = "BFGS", 
                                      ret_nll = TRUE, 
                                      server = TRUE,
                                      fit = FALSE,
                                      control = list(reltol = 1e-10, maxit = 100000)
)

# model_fit = optim(model_fit_lamE_betH_alpA$par,
#                   model_fit_lamE_betH_alpA$fn,
#                   model_fit_lamE_betH_alpA$gr,
#                   method = "BFGS",
#                   control = list(maxit = 1e5, reltol = 1e-9))
# std_errors = sdreport(model_fit_lamE_betH_alpA)

# Once again, this takes a while, but the object above can be optimized using the frequentist approach with the code above
stan_fit_lamE_betH_alpA = tmbstan(model_fit_lamE_betH_alpA, 
                                  cores = 5,
                                  init = "par", # initial values are all 0
                                  iter = 1250, # with the real dataset we increased this value until R-hat values were desirably low
                                  refresh = 5, 
                                  chains = 5)
