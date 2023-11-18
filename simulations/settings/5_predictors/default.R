# parameters_5_predictors.R
# Parameter definitions for a simulation with 5 predictors.

# Define the number of measurements, observations, basis functions, error standard deviation, and seed
predictors <- 5
true_predictors <- c(1,1,0,0,0)
measurements <- 100
observations <- 150
basis_functions = 5
intercept = 0
norder = 4
error_sd = 0.00
noise_sd = 0.00
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
  list(sig2 =0.5, rho = 0.5, decay_type = "matern"),
  list(sig2 =0.5, rho = 0.5, decay_type = "matern"),
  list(sig2 =0.5, rho = 0.5, decay_type = "matern"),
  list(sig2 =0.5, rho = 0.5, decay_type = "matern"),
  list(sig2 =0.5, rho = 0.5, decay_type = "matern")
)

beta_funcs <- list(
  function(t) 4 * exp(-abs(t - 1)),
  function(t) 2 * (1 - abs(t - 1)),
  function(t) 0 * t,
  function(t) 0 * t,
  function(t) 0 * t
)


time_domains = list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1),
  list(0, pi / 3),
  list(-2, 1)
)
