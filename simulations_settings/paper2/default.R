# Parameter definitions for a custom simulation

true_predictors <- c(1, 0)

predictors <- 2
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
  list(0, 300)

)

beta_funcs <- list(
  function(t) 5*dgamma(t/10, 3, 1/3),
  function(t) rep(0, length(t))
)

mu_funcs <- list(
  function(t) t,
  function(t) t
)

cov_funcs <- NULL

# # plot the gamma functions and give name
# t <- seq(0, 300, length.out=1000)
# plot(t, beta_funcs[[1]](t), type='l', col='blue', xlab='Time', ylab='Value', main=expression(paste("Beta_1: 5*", gamma, "(t/10, 3, 1/3)")))