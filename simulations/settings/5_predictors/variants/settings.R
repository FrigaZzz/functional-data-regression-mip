# parameters_5_predictors.R
# Parameter definitions for a simulation with 5 predictors.

# Define the number of measurements, observations, basis functions, error standard deviation, and seed
predictors <- 5
measurements <- 500
observations <- 500
basis_functions = 5
intercept = 0
norder = 4
error_sd = 0.05
seed = 1
# Define mu, covariance, and beta functions along with time domains
mu_funcs <- list(
  function(t) sin(2 * pi * t),
  function(t) 0.5 * (1 + sin(2 * pi * t)),
  function(t) 4 * t^2,
  function(t) 0.5 * (exp(-t) + exp(-2 * t) + exp(-3 * t) + exp(-4 * t) + exp(-5 * t)),
  function(t) 4/3 * log(1 + t^2)
)

cov_funcs <- list(
  list(sig2 = 0.2, rho = 0.11, decay_type = "matern"),
  list(sig2 = 0.22, rho = 0.21, decay_type = "matern"),
  list(sig2 = 0.42, rho = 0.31, decay_type = "matern"),
  list(sig2 = 0.2, rho = 0.11, decay_type = "matern"),
  list(sig2 = 0.22, rho = 0.21, decay_type = "matern")
)

beta_funcs <- list(
  function(t) 4 * exp(-abs(t - 1)),
  function(t) 2 * (1 - abs(t - 1)),
  function(t) 0 * t,
  function(t) 0 * t,
  function(t) 2 * t
)

time_domains <- list(
  seq(0, 1, length.out = measurements),
  seq(0, pi / 3, length.out = measurements), 
  seq(-1, 1, length.out = measurements), 
  seq(0, pi / 3, length.out = measurements), 
  seq(-2, 1, length.out = measurements)
)

