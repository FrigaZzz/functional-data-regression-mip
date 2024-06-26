# Parameter definitions for a custom simulation
# vector of affected betas by c value
c_val <- 0.0
true_predictors <- c(1, 1, 0, 1, 0, 0)

predictors <- 6
measurements <- 1000
observations <- 250
basis_functions <- 4
intercept <- 0
norder <- 4
noise_snr = c(100,100)
seed <- 1

simulation_type = "paper"
# cov_funcs is null only in the paper simulation
cov_funcs <- NULL



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
    rep(0, length(t)) 
  }
)






time_domains <- list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1),
  list(0, pi / 3),
  list(-2, 1),
  list(-1, 1)
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


# original
coef_specs <- list(
  '1' = list(
    a1 = list(type = "norm", mean = -4, sd = 3),
    a2 = list(type = "norm", mean = 7, sd = 1.5)
  ),
  '2' = list(
    b1 = list(type = "unif", min = 3, max = 7),
    b2 = list(type = "norm", mean = 0, sd = 1)
  ),
  '3' = list(
    c1 = list(type = "norm", mean = -3, sd = 1.2),
    c2 = list(type = "norm", mean = 2, sd = 0.5),
    c3 = list(type = "norm", mean = -2, sd = 1)
  ),
  '4' = list(
    d1 = list(type = "norm", mean = -2, sd = 1),
    d2 = list(type = "norm", mean = 3, sd = 1.5)
  ),
  '5' = list(
    e1 = list(type = "unif", min = 2, max = 7),
    e2 = list(type = "norm", mean = 2, sd = 0.4)
  ),
  '6' = list(
    f1 = list(type = "norm", mean = 4, sd = 2),
    f2 = list(type = "norm", mean = -3, sd = 0.5),
    f3 = list(type = "norm", mean = 1, sd = 1)
  )
)