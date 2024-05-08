source("FUNCTIONS.R")

require(TMB)
require(terra)
require(tidyverse)
require(rstan)
require(tmbstan)

SCALE_ELEV = 1000 # m to km
SCALE_DTWN = 1000 # m to km
SCALE_ALPHA_DIST = 1000 # m to km
SCALE_TSLV = 1 / 60 # change hours to minutes
N_SITES_TO_KEEP = 15 # In the real code this is 50, so it only includes the 50 closest sites to each camera (this is more than enough to capture all relevant variation according to the model's parameter estimates), but here we use 15 because there are only 15 sites anyways

SEASON_PARS = c(-5, 370, 0.01, 0.01) # implies no seasonal activity pattern

detections = read_csv("human_detections_SAMPLE.csv")
sites = read_csv("human_camera_effort_SAMPLE.csv")
all_site_data = read_csv("site_data_SAMPLE.csv")

## FILTER OUT SITES THAT ARE NOT IN STUDY AREA
all_site_data = all_site_data %>% dplyr::filter(in_study == 1)
detections = detections %>% dplyr::filter(site_name %in% all_site_data$site_name)
sites = sites %>% dplyr::filter(site_name %in% all_site_data$site_name)

# Add covariates for everything (scaling / squaring / etc)
detections = detections %>% mutate(intercept = 1,
                                   log_humans_per_day = log(humans_per_day + 1/3650),
                                   elevation_km = elevation / SCALE_ELEV, 
                                   slope2 = slope^2, cos_aspect = cos(aspect_radians), 
                                   sin_aspect = sin(aspect_radians), 
                                   dtown_km = dtown_m / SCALE_DTWN, 
                                   droad_km = droad_m / SCALE_DTWN, 
                                   dwater_km = dwater_m / SCALE_DTWN,
                                   snow_annual_2 = snow_annual ^ 2) %>% 
  mutate(elevation2_km = elevation_km ^ 2, elevation3_km = elevation_km ^ 3, 
         road01 = road01_paved + road01_unpaved, trail01 = trail01_formal + trail01_informal)

sites = sites %>% mutate(intercept = 1,
                         log_humans_per_day = log(humans_per_day + 1/3650),
                         elevation_km = elevation / SCALE_ELEV, 
                         slope2 = slope^2, 
                         cos_aspect = cos(aspect_radians), 
                         sin_aspect = sin(aspect_radians), 
                         dtown_km = dtown_m / SCALE_DTWN,
                         droad_km = droad_m / SCALE_DTWN, 
                         dwater_km = dwater_m / SCALE_DTWN,
                         snow_annual_2 = snow_annual ^ 2) %>%
  mutate(elevation2_km = elevation_km ^ 2, elevation3_km = elevation_km ^ 3, 
         road01 = road01_paved + road01_unpaved, trail01 = trail01_formal + trail01_informal)

site_histories = sites[, c("DateStart", "DateLastWorking")] %>% as.matrix
site_histories_full = sites %>% mutate(Site = site_name, 
                                       tstart = DateStart, 
                                       tend = DateLastWorking, 
                                       uptime = DateLastWorking - DateStart, .keep = "none")

all_site_data = all_site_data %>% mutate(log_humans_per_day = log(humans_per_day + 1/3650),
                                         trail01 = trail01_formal + trail01_informal)

# Which covariates are we using here? Notice I left out elevation and snow cover because there was trouble extrapolating from those variables
covariates_E = c("slope", "slope2", "cos_aspect", "sin_aspect", "barren", "herbaceous", "conif_open", "shrub", "droad_km", "dtown_km", "trail01_formal", "trail01_informal", "road01", "dwater_km")

# Prepare the necessary data
t_weights_log_sum = sum(log(seasonality_func(detections$t, SEASON_PARS, mod = TRUE)))
site_weights = int_seasonality_func(site_histories_full$tstart, site_histories_full$tend, SEASON_PARS, mod = TRUE)

t_covariates = cbind(1, as.matrix(detections[, covariates_E]))
t_covariates_sum = colSums(t_covariates)
site_covariates = cbind(1, as.matrix(sites[, covariates_E]))

covariates_H = c("humans_per_day")

# get distances between sites by converting to SpatVector
all_site_data_sp = vect(all_site_data, geom = c("easting", "northing"), crs = "+init=epsg:32611")
all_site_distances = distance(all_site_data_sp, all_site_data_sp) # Note that this includes the site itself (distance = 0) and it should

# Change distances between sites that are not either both on trail or both off trail to 200000, which is larger than any of the distances in the sites naturally, so they are not included in N_SITES_TO_KEEP
same_onoff_trail = outer(1:nrow(all_site_data), 1:nrow(all_site_data), FUN = function(r1, r2) {
  all_site_data$trail01[r1] == all_site_data$trail01[r2]
})
all_site_distances = ifelse(same_onoff_trail, all_site_distances, 200000)

connect_indices_by_site = t(apply(all_site_distances, 1, order))
connect_indices_by_site = connect_indices_by_site[, 1:N_SITES_TO_KEEP]

t_site_indices = sapply(1:nrow(detections), function(row) which(all_site_data$site_name == detections$site_name[row]))
site_site_indices = sapply(1:nrow(sites), function(row) which(all_site_data$site_name == sites$site_name[row]))

# these matrices are all integers representing the indices we need to draw from the distance and covariate objects
t_all_connections = connect_indices_by_site[t_site_indices, ]
site_all_connections = connect_indices_by_site[site_site_indices, ]

t_distances = sapply(1:nrow(detections), function(row) all_site_distances[t_site_indices[row], t_all_connections[row, ]]) %>% t
site_distances = sapply(1:nrow(sites), function(row) all_site_distances[site_site_indices[row], site_all_connections[row, ]]) %>% t

t_distances = t_distances / SCALE_ALPHA_DIST
site_distances = site_distances / SCALE_ALPHA_DIST

# get covariates (with the help of a little helper function)
get_cov_value = function(cov, connections = t_all_connections) {
  cov_vec = all_site_data[, cov, drop = TRUE]
  result = cov_vec[connections]
  dim(result) = dim(connections)
  result
}

t_spcov_values = sapply(covariates_H, FUN = get_cov_value, simplify = "array")
site_spcov_values = sapply(covariates_H, FUN = get_cov_value, simplify = "array", connections = site_all_connections)

# Comment this out since we are using log_humans_per_day
MIN_SPCOV_VALUE = 1e-7 # so the log doesn't go to 0 ever
t_spcov_values[t_spcov_values <= 0] = MIN_SPCOV_VALUE
site_spcov_values[site_spcov_values <= 0] = MIN_SPCOV_VALUE

# Remove the focal site (it's obviously not necessary for the regression in this case)
t_spcov_values = t_spcov_values[, -1, , drop = FALSE]
site_spcov_values = site_spcov_values[, -1, , drop = FALSE]
t_distances = t_distances[, -1]
site_distances = site_distances[, -1]

filename = "specialcovs_human_20231031"
compile(paste0(filename, ".cpp"))

dyn.load(dynlib(filename))

func_custom = MakeADFun(data = list(t_covariates = t_covariates,
                                    t_distances = t_distances,
                                    t_spcov_values = t_spcov_values[, , 1],
                                    site_weights = site_weights,
                                    site_covariates = site_covariates,
                                    site_distances = site_distances,
                                    site_spcov_values = site_spcov_values[, , 1],
                                    connection_ones = rep(1, N_SITES_TO_KEEP - 1)),
                        parameters = list(lambda_i = rep(0, ncol(t_covariates)),
                                          logalpha_i = 0,
                                          lognu_i = 0),
                        DLL = filename)

# model_fit = optim(func_custom$par,
#                   func_custom$fn,
#                   func_custom$gr,
#                   method = "BFGS",
#                   control = list(maxit = 1e5, reltol = 1e-9))
# std_errors = sdreport(func_custom)

# This takes a nontrivial amount of time, especially without access to parallel processing. Note that the above model object can also be optimized using the commented out lines above.
stan_fit = tmbstan(func_custom, 
                   cores = 5, # can turn this off if there isn't access to parallel processing
                   iter = 2500, # this is the # of iterations required to achieve sufficiently low R-hat values
                   init = 0, 
                   refresh = 5, 
                   chains = 5)
