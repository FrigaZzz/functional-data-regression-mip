# Parameter definitions for a custom simulation
# vector of affected betas by c value
c_val <- 0
true_predictors <- c(1, 1, 0, 1, 0, 0)

predictors <- 6
measurements <- 1000
observations <- 250
basis_functions <- 6
intercept <- 0
norder <- 4
noise_snr = c(100,100)
seed <- 1

simulation_type = "paper"

beta_funcs <- list(
  function(t) {
    sin(t)
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    if (c_val == 0) {
      rep(0, length(t))
    } else {
      -c_val * t^2
    }
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    if (c_val == 0) {
      rep(0, length(t))
    } else {
      c_val * sin(pi * t)
    }
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
  },
  function(t, args) {
    sin(2 * (t - args$d1)) + args$d2 * t
  },
  function(t, args) {
    args$e1 * cos(2 * t) + args$e2 * t
  },
  function(t, args) {
    args$f1 * exp(-t / 3) + args$f2 * t + args$f3
  }
)

# time_domains <- list(
#   list(0, 1),
#   list(0, 1),
#   list(0, 1),
#   list(0, 1),
#   list(0, 1),
#   list(0, 1)
# )

time_domains <- list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1),
  list(0, pi / 3),
  list(-2, 1),
  list(-1, 1)
)
# cov_funcs is null only in the paper simulation
cov_funcs <- NULL
