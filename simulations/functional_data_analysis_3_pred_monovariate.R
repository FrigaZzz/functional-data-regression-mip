
# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

source(here("src", "R", "data_simulator", "run_simulation.R"))

# Input params
params <- list(
  predictors = 3,
  observations = 500,
  measurements = 500,
  basis_functions = 5,
  intercept = 0,
  norder = 4,
  mu_funcs = list(
    function(t) sin(2 * pi * t),
    function(t) 4*t * (1 + cos(2 * pi * t)),
    function(t) 4 *t^2
  ),
  cov_funcs = list(
  # sig2, rho
    generate_covariance_function(0, 1,"matern"),
    generate_covariance_function(0, 1, "matern"),
    generate_covariance_function(0, 1,"matern")
  ),
  # t is the range of T values in which the function is evaluated
  beta_funcs <- list(
    function(t) cos(t),
    function(t)  sin(t),
    function(t) 0 * t # This assumes that the function returns a single scalar zero
  ),
  time_domain = seq(0, 1, length.out = measurements),

  error_sd =0,
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
true_predictors <- apply(B, 1, function(x) ifelse(all(x == 0), 0, 1))