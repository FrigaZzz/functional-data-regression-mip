
# Parameter definitions for a custom simulation
# vector of affected betas by c value
c_val <- 0
true_predictors <- c(1,1,1,1,1,0)

predictors <- 6
measurements = 50
observations = 100
basis_functions = 6
intercept = 0
norder = 4
error_sd = 0.05
noise_sd = 0.00
seed = 1
beta_funcs = list(
  function(t) {
    sin(t)
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    -c_val * t^2
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    c_val * sin(pi * t)
  },
  function(t) {
    0 *t # This function always returns 0 regardless of the input t
  }
)
coef_list <- list(
      list(a1 = rnorm(1, mean = -4, sd = 3), a2 = rnorm(1, mean = 7, sd = 1.5)),
      list(b1 = runif(1, min = 3, max = 7), b2 = rnorm(1, mean = 0, sd = 1)),
      list(c1 = rnorm(1, mean = -3, sd = sqrt(1.2^2)), c2 = rnorm(1, mean = 2, sd = sqrt(0.5^2)), c3 = rnorm(1, mean =    -2, sd = 1)),
      list(d1 = rnorm(1, mean = -2, sd = 1), d2 = rnorm(1, mean = 3, sd = sqrt(1.5^2))),
      list(e1 = runif(1, min = 2, max = 7), e2 = rnorm(1, mean = 2, sd = sqrt(0.4^2))),
      list(f1 = rnorm(1, mean = 4, sd = sqrt(2^2)), f2 = rnorm(1, mean = -3, sd = sqrt(0.5^2)), f3 = rnorm(1, mean = 1,     sd = 1))
  )
mu_funcs = list(
    function(t) {
        args <- coef_list[[1]] # Retrieve the fixed coefficients for this predictor
        cos(2 * pi * (t - args$a1)) + args$a2
    },
    function(t) {
      args <- coef_list[[2]]
        args$b1 * sin(pi * t) + args$b2
    },
    function(t) {
      args <- coef_list[[3]]
        args$c1 * t^3 + args$c2 * t^2 + args$c3 * t
    },
    function(t) {
      args <- coef_list[[4]]
        sin(2 * (t - args$d1)) + args$d2 * t
    },
    function(t) {
      args <- coef_list[[5]]
        args$e1 * cos(2 * t) + args$e2 * t
    },
    function(t) {
      args <- coef_list[[6]]
        args$f1 * exp(-t / 3) + args$f2 * t + args$f3
    }
)
time_domains = list(
  seq(0, 1, length.out = measurements),
  seq(0, pi / 3, length.out = measurements), 
  seq(-1, 1, length.out = measurements), 
  seq(0, pi / 3, length.out = measurements), 
  seq(-2, 1, length.out = measurements), 
  seq(-1, 1, length.out = measurements) 
)
# cov_funcs is null only in the paper simulation
cov_funcs = NULL

