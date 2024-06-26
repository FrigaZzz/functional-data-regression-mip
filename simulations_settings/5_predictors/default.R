# parameters_5_predictors.R
# Parameter definitions for a simulation with 5 predictors.
simulation_type = "cov"

# Define the number of measurements, observations, basis functions, error standard deviation, and seed
predictors <- 5
true_predictors <- c(1,1,0,0,1)
measurements <- 100
observations <- 150
basis_functions = 5
intercept = 0
norder = 4
noise_snr = c(100,100)
seed = 1
# Define mu, covariance, and beta functions along with time domains
mu_funcs <- list(
  function(t) sin(2 * pi * t),
  function(t) 0.5 * (1 + sin(2 * pi * t)),
  function(t) 4 * t^2,
  function(t) 0.5 * (exp(-t) + exp(-2 * t) + exp(-3 * t) + exp(-4 * t) + exp(-5 * t)),
  function(t) 1 /log(1 + t^2)
)

cov_funcs <- list(
  list(sig=0, rho = 0.1, decay_type = "exponential"),
  list(sig=0, rho = 0.1, decay_type = "matern"),
  list(sig=0, rho = 0.1, decay_type = "exponential"),
  list(sig=0, rho = 0.1, decay_type = "absolute_difference"),
  list(sig=0, rho = 0.1, decay_type = "absolute_difference")
)

beta_funcs <- list(
  function(t) 4 * exp(-abs(t - 1)),
  function(t) 2 * (1 - abs(t - 1)),
  function(t) 0 * t,
  function(t) 0 * t,
  function(t) 2 * t
)


time_domains = list(
  list(0, 3*pi),
  list(0, 5*pi / 3),
  list(-1,3* 1),
  list(0, 3*pi / 3),
  list(-2,3* 1)
)
