# Parameter definitions for a custom simulation
# vector of affected betas by c value
c_val <- 0
simulation_type = "paper"


true_predictors <- c(1, 1, 0, 0, 0, 0)

predictors <- 6
measurements <- 1000
observations <- 250
basis_functions <- 6
intercept <- 0
norder <- 4
noise_snr = c(100,100)
seed <- 1

coef_specs <- coef_specs_original 
beta_funcs <- list(
  function(t) {
    sin(t) + t/2 + cos(2 * pi * t)
  },
  function(t) {
    cos(t)
  },
  function(t) {
    rep(0, length(t))

  },
  function(t) {
    rep(0, length(t))

  },
  function(t) {
    rep(0, length(t))

  },
  function(t) {
    rep(0, length(t))
  }
)


mu_funcs <- list(
  function(t, args) {
    cos(2 * pi * (t - args$a1)) + args$a2
  },
  function(t, args) {
    0.5^(t-0.5) - cos(2 * pi * (t - args$b1)) + args$b2
  },
  function(t, args) {
    2 - cos(2 * pi * (t - args$c1)) + args$c2
  },
  function(t, args) {
    1/3 * cos(2 * pi * (t - args$d1)) + args$d2
  },
  function(t, args) {
    cos(3 * pi * (t - args$e1)) + args$e2
  },
  function(t, args) {
    cos(5 * pi * (t - args$f1)) + args$f2
  }
)

time_domains <- list(
  list(0, pi),
  list(0, pi),
  list(0, pi),
  list(0, pi),
  list(0, pi),
  list(0, pi)
)


# cov_funcs is null only in the paper simulation
cov_funcs <- NULL
