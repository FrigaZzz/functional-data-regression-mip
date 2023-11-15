

# run_simulation.R
# This script contains the function to run the simulation and returns all relevant outputs.

library(refund)
library(MASS)
library(fda)
library(here)

# Source the generic simulator script
source(here("src", "R", "generic_simulator", "simulate.R"))

# Define the function to run the simulation
run_simulation <- function(params) {
  # Ensure that the seed is set for reproducibility
  set.seed(params$seed)
  
  # Call the main analysis function with the parameters
  result <- do.call(run_functional_data_analysis, params)
  
  # Extract the required outputs
  output_list <- list(
    X = result$X,
    Y = result$Y,
    Z = result$Z,
    J = result$J,
    W = result$W,
    B = result$B,
    basis_objs = result$basis_objs,
    basis_values = result$basis_values,
    beta_point_values = result$beta_point_values,
    predictors = params$predictors,
    true_predictors = apply(result$B, 1, function(x) ifelse(all(x == 0), 0, 1))
  )
  
  # Return the output list
  return(output_list)
}

# Usage:
params <- list(
  predictors = predictors,
  measurements = measurements,
  observations = observations,
  basis_functions = basis_functions,
  intercept = intercept,
  norder = norder,
  error_sd = error_sd,
  seed = seed,
  mu_funcs = mu_funcs,
  cov_funcs = cov_funcs,
  beta_funcs = beta_funcs,
  time_domains = time_domains
)
outputs <- run_simulation(params)