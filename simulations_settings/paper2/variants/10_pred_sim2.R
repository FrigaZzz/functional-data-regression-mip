# Parameter definitions for a custom simulation

true_predictors <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)

predictors <- 10
measurements <- 1000
observations <- 250
basis_functions <- 4
intercept <- 0
norder <- 4
noise_snr <- c(100, 100)
seed <- 1

simulation_type <- "paper2"

time_domains <- list(
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300),
  list(0, 300)
)

beta_funcs <- list(
  function(t) (10 - (1.5 * 1)) * dgamma(t / 10, 3, 1 / 3),
  function(t) (10 - (1.5 * 2)) * dgamma(t / 10, 3, 1 / 3),
  function(t) (10 - (1.5 * 3)) * dgamma(t / 10, 3, 1 / 3),
  function(t) 0.6 * dexp(0.02 * t),
  function(t) 0.4 * dexp(0.01 * t),
  function(t) rep(0, length(t)),
  function(t) rep(0, length(t)),
  function(t) rep(0, length(t)),
  function(t) rep(0, length(t)),
  function(t) rep(0, length(t))
)

mu_funcs <- list(
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t,
  function(t) t
)

cov_funcs <- NULL

# # plot beta_funcs[3] and 4
# t <- seq(0, 300, length.out = 1000)
# # plot all of the betas functions in subplots
# par(mfrow = c(5, 1))
# plot(t, beta_funcs[[3]](t), type = "l", col = "blue", xlab = "Time", ylab = "Value", main = expression(paste("Beta_3: ", (10 - (1.5 * 3)), "*", gamma, "(t/10, 3, 1/3)")))
# plot(t, beta_funcs[[4]](t), type = "l", col = "blue", xlab = "Time", ylab = "Value", main = expression(paste("Beta_4: ", 0.6, "*", "exp(0.02 * t)")))
# plot(t, beta_funcs[[5]](t), type = "l", col = "blue", xlab = "Time", ylab = "Value", main = expression(paste("Beta_5: ", 0.4, "*", "exp(0.01 * t)")))

