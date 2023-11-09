# Load necessary libraries
library(refund)
library(MASS)
library(fda)
# Source utility files
source(here("src", "R", "data_simulator", "config.R")) # sets the utility paths

run_functional_data_analysis <- function(
    predictors, observations, measurements, basis_functions, intercept = 0, basis_measurements, norder,
    mu_funcs, cov_funcs, beta_funcs, error_sd, seed = 2000) {

  # Set seed for reproducibility
  set.seed(seed)

  # Define time domain
  time <- seq(0, 1, length.out = measurements)
  # Simulate functional features
  X_list <- simulate_functional_features(mu_funcs, cov_funcs, observations, time)

  # Define basis range and knots
  bm_range <- seq(0, 1, length.out = basis_measurements)
  knots <- bm_range[seq(1, basis_measurements, length.out = basis_functions - norder + 2)]

  # Create B-spline basis object
  basis_obj <- create.bspline.basis(rangeval = c(min(bm_range), max(bm_range)), nbasis = basis_functions, norder = norder, breaks = knots)

  # Expand functional data into B-spline basis
  expanded_data_list <- lapply(X_list, expand_functional_data, time = time, basis_obj = basis_obj)

  # Smooth beta coefficients into the B-spline basis
  expanded_beta_list <- lapply(beta_funcs, smooth_beta_function, time = time, basis_obj = basis_obj)

  dim(expanded_beta_list)
  # Initialize the array with correct dimensions for B coefficients
  B_array <- array(0, dim = c(length(expanded_beta_list), basis_functions))

  # Fill the array with the basis coefficients
  for (i in seq_along(expanded_beta_list)) {
    B_array[i, ] <- expanded_beta_list[[i]]
  }

  # Compute the J matrix (inner products for basis functions)
  J <- compute_J_matrix(basis_obj)

  # Initialize W_array with expanded data
  W_array <- array(unlist(expanded_data_list), dim = c(observations, predictors, basis_functions))

  # Compute Z matrix: Z = W %*% J (inner products for basis functions).
  # Dimensions: [observations, (1 +predictors) * basis_functions]
  Z_matrix <- compute_Z_matrix(W_array, J, predictors, basis_functions)

  Z_matrix <- array(Z_matrix, dim = c(observations,predictors,basis_functions))

  # Pre-generate all error terms
  all_errors <- rnorm(observations, mean = 0, sd = error_sd) # Generating error terms here

  # Compute Y values
  Y <- compute_Y_values(observations, predictors, basis_functions, B_array,intercept, Z_matrix, all_errors)

  # Return the computed Y values
  list(W = W_array, Z = Z_matrix, Y = Y, J = J, B = B_array, X = X_list, basis_obj = basis_obj, basis_values = eval.basis(bm_range, basis_obj))
}


