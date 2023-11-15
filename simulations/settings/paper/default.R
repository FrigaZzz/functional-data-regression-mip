
# Parameter definitions for a custom simulation
c = 0

predictors <- 6
measurements = 50
observations = 100
basis_functions = 6
intercept = 0
norder = 4
error_sd = 0.05
seed = 1
beta_funcs = list(
  function(t) {
    sin(t)
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    -c * t^2
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    c * sin(pi * t)
  },
  function(t) {
    0 *t # This function always returns 0 regardless of the input t
  }
)
 cov_funcs = list(
  # sig2, rho
    list(sig2 = 0.2, rho = 0.11,  decay_type = "matern"),
    list(sig2 = 0.22, rho = 0.21, decay_type = "matern"),
    list(sig2 = 0.42, rho = 0.31, decay_type ="matern"),
    list(sig2 = 0.12, rho = 0.11, decay_type = "matern"),
    list(sig2 = 0.32, rho = 0.01, decay_type ="matern"),
    list(sig2 = 0.32, rho = 0.01, decay_type ="matern")

  )
mu_funcs = list(
    function(t) {
        args <- list(a1 = rnorm(1, mean = -4,  sd = sqrt(3 ^ 2)), a2 = rnorm(1, mean = 7, sd = sqrt(1.5 ^ 2)))
        cos(2 * pi * (t - args$a1)) + args$a2
    },
    function(t) {
        args <- list(b1 = runif(1, min = 3, max = 7), b2 = rnorm(1, mean = 0, sd = 1))
        args$b1 * sin(pi * t) + args$b2
    },
    function(t) {
        args <- list(c1 = rnorm(1, mean = -3, sd = sqrt(1.2^2)), c2 = rnorm(1, mean = 2, sd = sqrt(0.5^2)), c3 = rnorm(1, mean = -2, sd = 1))
        args$c1 * t^3 + args$c2 * t^2 + args$c3 * t
    },
    function(t) {
        args <- list(d1 = rnorm(1, mean = -2, sd = 1), d2 = rnorm(1, mean = 3, sd = sqrt(1.5^2)))
        sin(2 * (t - args$d1)) + args$d2 * t
    },
    function(t) {
        args <- list(e1 = runif(1, min = 2, max = 7), e2 = rnorm(1, mean = 2, sd = sqrt(0.4^2)))
        args$e1 * cos(2 * t) + args$e2 * t
    },
    function(t) {
        args <- list(f1 = rnorm(1, mean = 4, sd = sqrt(2^2)), f2 = rnorm(1, mean = -3, sd = sqrt(0.5^2)), f3 = rnorm(1, mean = 1, sd = 1))
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

