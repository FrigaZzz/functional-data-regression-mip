# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)
# Source utility files
source(here("src", "R", "data_simulator", "config.R")) # sets the utility paths

run_functional_data_analysis <- function(
    predictors, observations, measurements, basis_functions, intercept = 0, norder,
    mu_funcs, cov_funcs, beta_funcs, time_domain, error_sd, seed = 2000) {

  # Set seed for reproducibility
  set.seed(seed)

  # Define time domain
  # 1. Simulate functional features
  X <- simulate_functional_features(mu_funcs, cov_funcs, observations, time_domain)
 

  # 2. Create B-spline basis object
  basis_obj <- create_1_basis(basis_functions, time_domain, norder)

  # 3. Smooth beta coefficients into the B-spline basis
  result <- smooth_betas(beta_funcs,basis_functions, time_domain, basis_obj)
  B_array <- result$beta_matrix
  basis_values <- result$basis_values
  beta_point_values <- result$beta_point_values


  # 4. Compute the J matrix (inner products for basis functions)
  J <- compute_J_matrix(basis_obj)

  # 5. W_array: Expand X functional data into B-spline basis
  # Initialize W_array with expanded data
  W_array <- compute_W_matrix(X, basis_functions, time_domain, basis_obj)

  # 6. Compute Z matrix: Z = W %*% J 
  Z_matrix <- compute_Z_matrix(W_array, J, predictors, basis_functions)

  # 7. Compute Y values
  Y <- compute_Y_values(observations, predictors, basis_functions, B_array,intercept, Z_matrix)

  # 8. Apply error terms to Y values
  noise = compute_noise(Y, error_sd)
  Y <- Y + noise

  # Return the computed Y values
  list(W = W_array, Z = Z_matrix, Y = Y, J = J, B = B_array, X = X, basis_obj = basis_obj, basis_values = basis_values,beta_point_values = beta_point_values)
}


