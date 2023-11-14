
# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

source(here("src", "R", "paper_simulator", "paper_sim.R"))

c = 0
# measurements = 100
# observations = 500
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


# Input params
params <- list(
  predictors = 6,
  observations = observations,
  measurements= measurements,
  basis_functions = 6,
  intercept = 0,
  norder = 4,
  # t is the range of T values in which the function is evaluated
  mu_funcs= mu_funcs,
  beta_funcs= beta_funcs,
  time_domains= time_domains,
  seed = 1
)
# Call the main analysis function with the parameters
result <- do.call(run_functional_data_analysis, params)
# Print the first few rows of Y
# unpack the list of outputs 
X <- result$X
Y <- result$Y
Z <- result$Z
J <- result$J
W <- result$W
B <- result$B
basis_obj <- result$basis_obj
basis_values <- result$basis_values
beta_point_values <- result$beta_point_values
predictors <- params$predictors
true_predictors <- apply(B, 1, function(x) ifelse(all(x == 0), 0, 1))