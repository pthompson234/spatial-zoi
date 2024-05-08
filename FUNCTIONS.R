# Fit a intensity-detection model to camera trap data
#
# data: a list containing necessary information for model fitting
# init_par: either NULL or a numeric vector of the appropriate length, representing a refined starting position for the optimizer. Typically only useful when called by other functions (e.g., fit_ID_MonteCarloIntegral)
# sp_logit: logical; only matters if running "specialcovs" model; determines which version of the model is run
# fit: logical; do we fit the model or just return the function?
# server: logical; are we using this on a remote directory or the local directory
# bayesian: logical; do we fit the model using Stan or MLE?
# ret_nll: logical; do we return the TMB Object containing the likelihood function (and other things)?
# custom_filename: character: do we want to do TMB with a different file than usual?
# TMB_compile_framework: see ?TMB::compile() for more info
# internal_AD: logical; equivalent to "intern" argument for TMB::MakeADFun(), usually not necessary so set to FALSE by default
# ...: additional arguments, either to optimizer or to Stan sampler
#
# Returns a data.frame containing model fitting information (or if ret_nll is TRUE, a list with a data.frame and the TMB object)
fit_ID_TMB = function(data, 
                      init_par = NULL, 
                      sp_logit = FALSE, 
                      fit = TRUE,
                      server = FALSE,
                      bayesian = FALSE, 
                      ret_nll = FALSE,
                      custom_filename = "",
                      TMB_compile_framework = getOption("tmb.ad.framework"),
                      internal_AD = FALSE,
                      ...) {
  
  require(TMB)
  
  if (!server) {
    old_wd = getwd()
    on.exit(setwd(old_wd))
    setwd("~/School/Postdocs/CanmoreCorridor_20221017/code/cameras")
  }
  
  tmb_var_names = names(data)
  random = NULL # set to this unless we want to re-assign it inside the if statement
  
  if ("t_alcov_values" %in% tmb_var_names) {
    # running the special_covs model where alpha's also depend on covariates
    filename = "specialcovs_20231016"
    
    n_lambdas = length(data$t_covariates_sum)
    n_betas = dim(data$t_spcov_values)[3]
    n_beta_covs = ncol(data$t_alcov_values)
    
    parameters = list(lambda_i = numeric(n_lambdas), beta_i = matrix(0, n_betas, n_beta_covs), logalpha_i = matrix(0, n_betas, n_beta_covs))
  } else {
    stop("Unsupported model type (for now)")
  }
  
  if (nchar(custom_filename) > 0) filename = custom_filename
  
  # for debugging purposes
  # message("Current working directory: ", getwd())
  compile(paste0(filename, ".cpp"), framework = TMB_compile_framework)
  dyn.load(dynlib(filename))
  
  func = MakeADFun(data = data, parameters = parameters, random = random, DLL = filename, intern = internal_AD)
  
  if (!fit) return(func)
  if (is.null(init_par)) init_par = func$par
  
  # for now leave this in
  # if (bayesian) stop("Bayesian inference currently not implemented. Sorry; working on it!")
  if (bayesian) {
    require(tmbstan)
    return(tmbstan(func, ...))
  }
  
  # this_data <<- data
  # func <<- func
  # stop("Interrupted on purpose :)")
  
  model_fit = optim(init_par, func$fn, func$gr, ...)
  
  std_errors = sdreport(func)
  convergence_message = "Convergence error"
  if (model_fit$convergence == 1) convergence_message = "Iteration limit reached"
  if (model_fit$convergence == 0) convergence_message = "Converged"
  
  # for debugging purposes
  # return(list(model = model_fit, sdreport = std_errors))
  
  par_names = rownames(std_errors$cov.fixed)
  if (is.null(par_names)) par_names = paste0("Par ", 1:nrow(std_errors$cov.fixed))
  
  result = data.frame(Par = par_names, estimate = model_fit$par, se = std_errors$sd, NLL = model_fit$value, msg = convergence_message)
  res_df = result %>% mutate(CL = estimate - 1.96 * se, CU = estimate + 1.96*se)
  
  if (ret_nll) return(list(model = res_df, func = func))
  list(model = res_df)
  
}

# Calculate the seasonality weights (the difference between two appropriately shifted tanh() terms)
#
# t: numeric vector; times at which the weight should be calculated
# pars: numeric vector of length 4; parameters (theta1, theta2, phi1, phi2) for the model. Note that if phi1 and phi2 both equal 0 then it changes to a step function
# mod: logical; do we mod all t by 365 (i.e., if we are given times from different years?)
#
# Returns a numeric vector with seasonality values (which should range from 0 to 1)
seasonality_func = function(t, pars, mod = FALSE) {
  
  if (mod) t = t %% 365
  
  0.5 * (tanh((t - pars[1]) / pars[3]) - tanh((t - pars[2]) / pars[4]))
}

# Calculate the integral of our seasonality function (the difference between two appropriately shifted tanh() terms)
#
# t1: numeric vector; lower bounds for integration
# t2: numeric vector; upper bounds for integration
# pars: numeric vector of length 4; parameters (theta1, theta2, phi1, phi2) for the model. Note that if phi1 and phi2 both equal 0 then it changes to a step function
# mod: logical; passes to seasonality_func
#
# Returns a numeric vector with seasonality values, integrated (which should all be greater than 0)
int_seasonality_func = function(t1, t2, pars, mod = FALSE) {
  
  # While seasonality_func does have an analytical antiderivative, it involves lots of exponentials and can quickly become Inf or NaN. The numerical integral isn't super fast but for the purposes with which we are using this function, that's OK.
  Vectorize(function(x1, x2) {
    integrate(seasonality_func, x1, x2, pars = pars, mod = mod)$value
  }, vectorize.args = c("x1", "x2"))(t1, t2)
  
}