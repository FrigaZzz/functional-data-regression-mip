
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
  basis_measurements = 300,
  norder = 4,
  mu_funcs = list(
    function(t) sin(2 * pi * t),
    function(t) 0.5 * (1 + sin(2 * pi * t)),
    function(t) 40 *t^2
  ),
  cov_funcs = list(
  # sig2, rho
    generate_covariance_function(1, 0.1,"exponential"),
    generate_covariance_function(0.8, 0.3, "matern"),
    generate_covariance_function(0.5, 0.2,"exponential")
  ),
  # t is the range of T values in which the function is evaluated
  beta_funcs <- list(
    function(t) 4 * exp(-abs(t - 1)),
    function(t) 2 * (1 - abs(t - 1)),
    function(t) 0 * t # This assumes that the function returns a single scalar zero
  ),
  error_sd = 0.5,
  seed = 2000
)

# Call the main analysis function with the parameters
result <- do.call(run_functional_data_analysis, params)
# Print the first few rows of Y
# unpack the list of outputs 
Y <- result$Y
Z <- result$Z
J <- result$J
W <- result$W
B <- result$B
basis_obj <- result$basis_obj
basis_values <- result$basis_values