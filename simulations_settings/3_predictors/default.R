# parameters_3_predictors.R
simulation_type = "cov"

# Parameter definitions for a simulation with 3 predictors.

beta_funcs <- list(
  function(t) cos(t),
  function(t) sin(t),
  function(t) rep(0, length(t))
)

mu_funcs <- list(
  function(t) sin(pi * t),
  function(t) 2*t,
  function(t) cos(pi * t)
)

time_domains <- list(
  list(0, 2*pi),
  list(0, 2*pi),
  list(0, 2*pi)
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

cov_funcs <- list(
  list(sig2 = 0.5, rho = 1, decay_type = "matern"),
  list(sig2 = 0.5, rho = 1, decay_type =  "matern"),
  list(sig2 = 0.5, rho = 1, decay_type =  "matern")
)


# # plot each beta_funcs 
# # Create sumplot 3 plot based on the beta_funcs and time_domains.
# t <- seq(0, 2 * pi, length.out = 100)
# par(mfrow = c(1, 3))
# for (i in 1:3) {
#   plot(t, beta_funcs[[i]](t), type = "l", main = paste("beta", i), xlab = "t", ylab = "beta")
#   abline(v = time_domains[[i]], col = "red")
# }
# par(mfrow = c(1, 1))

# Create sumplot 3 plot based on the mu_funcs and time_domains. Only show them in corresponding time domain.
# t3 <- seq(0, 2*pi,length.out = 100)
# t1 <- seq(0, 2*pi,length.out = 100)
# t2 <- seq(0, 2*pi,, length.out = 100)
# par(mfrow = c(1, 3))
# # plo the beta functions
# plot(t1, beta_funcs[[1]](t1), type = "l", main = "b1", xlab = "t", ylab = "b")
# plot(t2, beta_funcs[[2]](t2), type = "l", main = "b2", xlab = "t", ylab = "b")
# plot(t3, beta_funcs[[3]](t3), type = "l", main = "b3", xlab = "t", ylab = "b")
# par(mfrow = c(1, 1))






