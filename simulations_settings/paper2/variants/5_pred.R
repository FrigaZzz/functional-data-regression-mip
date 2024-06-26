# Parameter definitions for a custom simulation

true_predictors <- c(1, 1, 1, 1,0)

predictors <- 5
measurements <- 1000
observations <- 250
basis_functions <- 4
intercept <- 0
norder <- 4
noise_snr = c(100,100)
seed <- 1

simulation_type = "paper2"

time_domains <- list(
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300)


)

beta_funcs <- list(
  function(t) 5*dgamma(t/10, 3, 1/3),
  function(t) 2*dgamma(t/10, 3, 2/3),
  function(t) 1*dgamma(2*t/10, 3, 1/3),
  function(t) dgamma(t/20, 3, 1/3),
  function(t) rep(0, length(t))
)

mu_funcs <- list(
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t
)

cov_funcs <- NULL
