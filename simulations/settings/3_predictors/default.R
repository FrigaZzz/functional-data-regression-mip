# parameters_3_predictors.R
# Parameter definitions for a simulation with 3 predictors.

beta_funcs <- list(
  function(t) cos(t),
  function(t) sin(t),
  function(t) 0 * t
)
predictors <- 3
measurements <- 500
observations <- 500
basis_functions = 5
intercept = 0
norder = 4
error_sd = 0.05
noise_sd = 0.00
seed = 1
mu_funcs <- list(
  function(t) sin(2 * pi * t),
  function(t) 4 * t * (1 + cos(2 * pi * t)),
  function(t) 4 * t^2
)

cov_funcs <- list(
  list(sig2 = 0.2, rho = 0.11, decay_type = "matern"),
  list(sig2 = 0.22, rho = 0.21, decay_type = "matern"),
  list(sig2 = 0.42, rho = 0.31, decay_type = "matern")
)
time_domains <- list(
  seq(0, 1, length.out = measurements),
  seq(0, pi / 3, length.out = measurements), 
  seq(-1, 1, length.out = measurements)
)

