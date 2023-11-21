

# run_simulation.R
# This script contains the function to run the simulation and returns all relevant outputs.

library(refund)
library(MASS)
library(fda)
library(here)

# Source the generic simulator script
source(here("src", "R", "generic_simulator", "simulate_main.R"))
source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))


# Define the function to run the simulation
run_simulation <- function(params) {
  # Ensure that the seed is set for reproducibility
  
  # Call the main analysis function with the parameters
  result <- do.call(generate_data, params[c("predictors", "observations", "measurements", "basis_functions", "intercept", "norder", "mu_funcs", "beta_funcs","time_domains", "cov_funcs", "seed","noise_snr","simulation_type")] )
    
  # Return the output list
  return(result)
}


simulation_name = "3_predictors"
simulation_settings_file = "default"
# Required inputs before running the simulation!!!
inputs  <- load_simulation_settings(simulation_name, simulation_settings_file)
param_name = "measurements"
inputs[[param_name]] <- 500
inputs$observations <- 120
inputs$basis_functions <- 6
time_domains_eval <- lapply(inputs$time_domains, function(domain) {
           seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
       })
inputs$time_domains <- time_domains_eval
inputs$noise_snr <- c(FALSE,100)
inputs$simulation_type <- "non paper"
outputs <- run_simulation(inputs);

