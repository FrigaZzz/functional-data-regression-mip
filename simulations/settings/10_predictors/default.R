# parameters_10_predictors.R
# This file defines parameters for a simulation with 10 predictors.
true_predictors <- c(1,1,0,0,0,0,0,0,0,1)

# Define the number of measurements and observations
predictors <- 10
measurements <- 500
observations <- 500
basis_functions = 5
intercept = 0
noise_snr = c(100,100)
norder = 4
seed = 1
mu_funcs = list(
    function(t) sin(2 * pi * t),
    function(t) 0.5 * (1 + sin(2 * pi * t)),
    function(t) 4 *t^2,
    function(t) 3 * t^3 - t^2 + 2,
    function(t) tanh(5 * t),
    function(t) exp(-(t^2)/2),
    function(t) abs(t) + t,
    function(t) 2 * (t - floor(t + 0.5)),
    function(t) 1 / (1 + t^2), 
    function(t) sin(0.5 * sin(pi * t)) * exp(0.2 * t)
  )

  cov_funcs = list(
  list(sig2 = 0.2, rho = 0.11, decay_type = "matern"),
  list(sig2 = 0.22, rho = 0.21, decay_type = "matern"),
  list(sig2 = 0.42, rho = 0.31, decay_type = "matern"),
  list(sig2 = 0.2, rho = 0.11, decay_type = "matern"),
  list(sig2 = 0.22, rho = 0.21, decay_type = "matern"),
  list(sig2 = 0.2, rho = 0.11, decay_type = "matern"),
  list(sig2 = 0.22, rho = 0.21, decay_type = "matern"),
  list(sig2 = 0.42, rho = 0.31, decay_type = "matern"),
  list(sig2 = 0.2, rho = 0.11, decay_type = "matern"),
  list(sig2 = 0.22, rho = 0.21, decay_type = "matern"),
  )
  beta_funcs = list(
    function(t) 4 * exp(-abs(t - 1)),
    function(t) 2 * (1 - abs(t - 1)),
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 4^t + t  
  )

  time_domains = list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1),
  list(0, pi / 3),
  list(-2, 1),
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1),
  list(0, pi / 3),
  list(-2, 1)
)

