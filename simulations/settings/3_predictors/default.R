# parameters_3_predictors.R
# Parameter definitions for a simulation with 3 predictors.

beta_funcs <- list(
  function(t) cos(t),
  function(t) sin(t),
  function(t) 0 * t
)
true_predictors <- c(1,1,0)

predictors <- 3
measurements <- 200
observations <- 100
basis_functions = 5
intercept = 0
norder = 4
noise_snr = c(100,100)
seed = 1
mu_funcs <- list(
  function(t) sin(2 * pi * t),
  function(t) 4 * t * (1 + cos(2 * pi * t)),
  function(t) 4 * t^2
)

cov_funcs <- list(
  list(sig2 = 0.5, rho = 1, decay_type = "matern"),
  list(sig2 = 0.5, rho = 1, decay_type =  "matern"),
  list(sig2 = 0.5, rho = 1, decay_type =  "matern")
)

time_domains = list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1)
)
