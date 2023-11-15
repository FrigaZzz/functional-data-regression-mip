

# run_simulation.R
# This script contains the function to run the simulation and returns all relevant outputs.

library(refund)
library(MASS)
library(fda)
library(here)

# Source the generic simulator script
source(here("src", "R", "generic_simulator", "simulate.R"))
source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))


# Define the function to run the simulation
run_simulation <- function(params) {
  # Ensure that the seed is set for reproducibility
  set.seed(params$seed)
  
  # Call the main analysis function with the parameters
  result <- do.call(generate_data, params[c("predictors", "observations", "measurements", "basis_functions", "intercept", "norder", "mu_funcs", "beta_funcs","time_domains", "cov_funcs", "error_sd", "seed","noise_sd")] )
  
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
    beta_point_values = result$beta_point_values
  )
  
  # Return the output list
  return(output_list)
}


# simulation_name = "3_predictors"
# simulation_settings_file = "1_active"
# Required inputs before running the simulation!!!
inputs  <- load_simulation_settings(simulation_name, simulation_settings_file)
