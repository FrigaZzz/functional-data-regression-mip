# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

# Source utility files
source(here("src","R", "generic_simulator", "config.R")) # sets the utility path
source(here("src", "R",  "generic_simulator",   "simulation" , "robust.R"))

# Unified function
generate_data <- function(
    predictors, observations, measurements, basis_functions, intercept = 0, norder, noise_snr, 
    mu_funcs= NULL, beta_funcs = NULL, time_domains = NULL, cov_funcs = NULL, seed = 2000,
    simulation_type = "not paper") {

  set.seed(seed)
  # debug print all input parameters
  print(paste("predictors:", predictors))
  print(paste("observations:", observations))
  print(paste("measurements:", measurements))
  print(paste("basis_functions:", basis_functions))
  print(paste("intercept:", intercept))
  

  # Call the appropriate data generation function based on simulation type
  data = simulate_data_robust(predictors, observations, measurements) 

  # Extract U, X, Y from the returned data
  U <- data$U
  X <- data$X
  Y <- data$Y
  # create beta values from each time domain
  beta_point_values <- data$beta_values
  time_domains <- rep(list(seq(0, 1, length.out = measurements)), predictors) # grid points of Xs

  # # Remaining processing steps
  basis_objs <- create_basis(basis_functions, time_domains, norder, predictors)
  result <- smooth_betas_generic(beta_point_values, basis_functions, time_domains, basis_objs)
  Beta_matrix <- result$beta_matrix
  basis_values <- result$basis_values
  J <- compute_J_matrix_generic(basis_objs, predictors, basis_functions)
  W <- compute_W_matrix_generic(X, basis_functions, time_domains, basis_objs)
  Z_matrix <- compute_Z_matrix_generic(W, J, predictors, basis_functions)

  list(W = W, Z = Z_matrix, Y = Y, J = J, B = Beta_matrix, U = U, X = X, basis_objs = basis_objs, basis_values = basis_values, beta_point_values = beta_point_values)
}
