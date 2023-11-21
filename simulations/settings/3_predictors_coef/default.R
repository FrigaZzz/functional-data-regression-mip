# Parameter definitions for a custom simulation
# vector of affected betas by c value
c_val <- 0
true_predictors <- c(1, 1, 0)

predictors <- 3
measurements <- 50
observations <- 100
basis_functions <- 6
intercept <- 0
norder <- 4
noise_snr = c(TRUE,TRUE)
seed <- 1

simulation_type = "paper"

beta_funcs <- list(
  function(t) {
    sin(2 * t * pi)
  },
  function(t) {
    sin(pi * t)
  },
  function(t) {
    rep(0, length(t)) # This function always returns 0 regardless of the input t
  }
)


mu_funcs <- list(
  function(t, args) {
    cos(2 * pi * (t - args$a1)) + args$a2
  },
  function(t, args) {
    args$b1 * sin(pi * t) + args$b2
  },
  function(t, args) {
    args$c1 * t^3 + args$c2 * t^2 + args$c3 * t
  }
)


time_domains <- list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1)
)
cov_funcs <- NULL
