# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

# Source utility files
source(here("src","R", "generic_simulator", "config.R")) # sets the utility path

# Unified function
generate_data <- function(
    predictors, observations, measurements, basis_functions, intercept = 0, norder, noise_snr, 
    mu_funcs, beta_funcs, time_domains, cov_funcs = NULL, seed = 2000,
    simulation_type = "not paper") {

  # set.seed(seed)

  # Call the appropriate data generation function based on simulation type
  if (simulation_type == "paper") {
    data <- simulate_paper_data(mu_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr)
  } else{
    data <- simulate_cov_data(mu_funcs, cov_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr) 
  }

  # Extract U, X, Y from the returned data
  U <- data$U
  X <- data$X
  Y <- data$Y

  # Remaining processing steps
  basis_objs <- create_basis(basis_functions, time_domains, norder, predictors)
  result <- smooth_betas_generic(beta_funcs, basis_functions, time_domains, basis_objs)
  Beta_matrix <- result$beta_matrix
  basis_values <- result$basis_values
  beta_point_values <- result$beta_point_values
  J <- compute_J_matrix_generic(basis_objs, predictors, basis_functions)
  W <- compute_W_matrix_generic(X, basis_functions, time_domains, basis_objs)
  Z_matrix <- compute_Z_matrix_generic(W, J, predictors, basis_functions)

  list(W = W, Z = Z_matrix, Y = Y, J = J, B = Beta_matrix, U = U, X = X, basis_objs = basis_objs, basis_values = basis_values, beta_point_values = beta_point_values)
}
