
# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

source(here("src", "R", "data_simulator", "run_simulation.R"))
observations <- 500
measurements <- 500
basis_functions <- 5
error_sd <-0.001
seed  <- 1


# Input params
params <- list(
  predictors = 5,  # You might want to make this dynamic as well
  observations = observations,
  measurements = measurements,
  basis_functions = basis_functions,
  intercept = 0,
  norder = 4,
  mu_funcs = list(
    function(t) sin(2 * pi * t),
    function(t) 0.5 * (1 + sin(2 * pi * t)),
    function(t) 4 *t^2,
    function(t) 0.5 * (exp(-t) + exp(-2 * t) + exp(-3 * t) + exp(-4 * t) + exp(-5 * t)),
    function(t) 4/3 * log(1 + t^2)
  ),
  cov_funcs = list(
  # sig2, rho
    generate_covariance_function(1, 1,"matern"),
    generate_covariance_function(1, 1, "matern"),
    generate_covariance_function(1, 1,"matern"),
    generate_covariance_function(1, 1, "matern"),
    generate_covariance_function(1, 1,"matern")
  ),
  # t is the range of T values in which the function is evaluated
  beta_funcs <- list(
    function(t) 4 * exp(-abs(t - 1)),
    function(t) 2 * (1 - abs(t - 1)),
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 2 * t # This assumes that the function returns a single scalar zero
  ),
  error_sd = error_sd,
  seed = seed
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
predictors <- params$predictors
# compute true_predictors as a vector of 0s and 1s based on this: if a row of B is all 0s, then the corresponding predictor is not a true predictor
# this is done for each row
true_predictors <- apply(B, 1, function(x) ifelse(all(x == 0), 0, 1))


