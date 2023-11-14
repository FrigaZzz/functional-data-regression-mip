
# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

source(here("src", "R", "data_simulator", "run_simulation.R"))

# Input params
params <- list(
  predictors = 10,
  observations = 500,
  measurements = 500,
  basis_functions = 5,
  intercept = 0,
  basis_measurements = 300,
  norder = 4,
  mu_funcs = list(
    function(t) sin(2 * pi * t),
    function(t) 0.5 * (1 + sin(2 * pi * t)),
    function(t) 4 *t^2,
    function(t) 3 * t^3 - t^2 + 2,
    function(t) tanh(5 * t),
    function(t) exp(-(t^2)/2),
    function(t) abs(t) + t,
    function(t) 2 * (t - floor(t + 0.5)),
    function(t) 1 / (1 + t^2), 
    function(t) sin(0.5 * sin(pi * t)) * exp(0.2 * t)
  ),
  cov_funcs = list(
    generate_covariance_function(1, 0.1),
    generate_covariance_function(0.8, 0.3, "matern"),
    generate_covariance_function(0.3, 0.2),
    generate_covariance_function(1, 0.3),
    generate_covariance_function(0.3, 0.3, "matern"),
    generate_covariance_function(0.3, 0.3),
    generate_covariance_function(0.2, 0.3),
    generate_covariance_function(0.8, 0.2),
    generate_covariance_function(0.6, 0.2),
    generate_covariance_function(1, 0.1)
  ),
  beta_funcs = list(
    function(t) 4 * exp(-abs(t - 1)),
    function(t) 2 * (1 - abs(t - 1)),
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t # This assumes that the function returns a single scalar zero
  ),
  time_domain = seq(0, 1, length.out = measurements),
  #   time_domains <- replicate(num_elements, seq(0, 1, length.out = measurements), simplify = FALSE),
  error_sd = 0.05,
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
beta_point_values <- result$beta_point_values