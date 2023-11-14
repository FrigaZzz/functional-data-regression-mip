# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)
# Source utility files
source(here("src","R", "generic_simulator", "config.R")) # sets the utility path


run_functional_data_analysis <- function(
    predictors, observations, measurements , basis_functions, intercept = 0, norder,
    mu_funcs, beta_funcs, time_domains, seed = 2000, noise_sd = 0) {

  # Set seed for reproducibility
  set.seed(seed)


  # 1. Simulate functional features
  X <- simulate_data_paper(observations, predictors, time_domains, mu_funcs)
  # Apply amplitude normalization

  # X <- apply_amplitude_norm(X, observations, predictors)
  
  # 2. Create B-spline basis object
  basis_objs <- create_basis(basis_functions, time_domains, norder, predictors)

  # 3. Smooth beta coefficients into the B-spline basis
  result <- smooth_betas_paper(beta_funcs,basis_functions, time_domains, basis_objs)
  Beta_matrix <- result$beta_matrix
  basis_values <- result$basis_values
  beta_point_values <- result$beta_point_values

  # 4. Compute the J matrix (inner products for basis functions)
  J <- compute_J_matrix_paper(basis_objs, predictors, basis_functions)

  # 5. W_array: Expand X functional data into B-spline basis
  W <- compute_W_matrix_paper(X, basis_functions, time_domains, basis_objs)

  # 6. Compute Z matrix: Z = W %*% J 
  Z_matrix <- compute_Z_matrix_paper(W, J, predictors, basis_functions)

  # 7. Compute Y values
  Y <- compute_Y_values(observations, predictors, basis_functions, Beta_matrix, intercept, Z_matrix)

  # 8. Apply error terms to Y values
  
  # Y <- Y + compute_amplitude_norm(Y)
  # Y <- Y + compute_noise(Y,noise_sd)


  # Return the computed Y values
  list(W = W, Z = Z_matrix, Y = Y, J = J, B = Beta_matrix, X = X, basis_obj = basis_objs, basis_values = basis_values,beta_point_values = beta_point_values)
}


