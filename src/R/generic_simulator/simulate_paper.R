
# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

# Source utility files
source(here("src","R", "generic_simulator", "config.R")) # sets the utility path

#' 
#' This function simulates functional features.
#' It takes in predictors, observations, measurements, basis functions, intercept, norder, mu_funcs, beta_funcs,
#' time_domains, cov_funcs, error_sd, seed, and noise_sd as inputs and returns a list of computed Y values.
#'
#' @param predictors A numeric vector of predictor values
#' @param observations A numeric vector of observation values
#' @param measurements A numeric vector of measurement values
#' @param basis_functions A character vector of basis functions
#' @param intercept A numeric value for the intercept
#' @param norder An integer value for the order of the B-spline basis
#' @param mu_funcs A list of functions for the mean of the functional features
#' @param beta_funcs A list of functions for the beta coefficients
#' @param time_domains A list of time domains for the functional features
#' @param cov_funcs A list of covariance functions for the functional features
#' @param seed An integer value for the random seed
#' @param noise_snr Defines the amount of noise to add to the data 
#'
#' @return A list containing the computed Y values

generate_data <- function(
    predictors, observations, measurements , basis_functions, intercept = 0, norder,
    mu_funcs, beta_funcs, time_domains, cov_funcs=NULL,  seed = 2000, noise_snr = c(100,100)) {

  # Set seed for reproducibility
  set.seed(seed)
  
  # 1. Simulate input data
  U <- NULL
  X <- NULL

  # Simulate true predictors
  U <- simulate_true_predictors_Ut(mu_funcs, observations, time_domains) 
  Betas <- create_beta_curves(beta_funcs,  time_domains)

  # Apply amplitude normalization to get observed features
  # X <-  simulate_observations_Xt(U) # adds the modulation
  # X <-  apply_snr_to_X(X, noise_snr[1])
  # copy U in X
  X <- U
  Y <- compute_Y_values(U, Betas, observations, predictors,time_domains, intercept)$Y

  # Apply error terms to Y values
  # Y <- Y + compute_amplitude_norm(Y) # adds the modulation


  # 2. Create B-spline basis object
  basis_objs <- create_basis(basis_functions, time_domains, norder, predictors)
  
  # 3. Smooth beta coefficients into the B-spline basis
  result <- smooth_betas_generic(beta_funcs,basis_functions, time_domains, basis_objs)
  Beta_matrix <- result$beta_matrix
  basis_values <- result$basis_values
  beta_point_values <- result$beta_point_values

  # 4. Compute the J matrix (inner products for basis functions)
  J <- compute_J_matrix_generic(basis_objs, predictors, basis_functions)

  # 5. W_array: Expand X functional data into B-spline basis
  W <- compute_W_matrix_generic(U, basis_functions, time_domains, basis_objs)

  # 6. Compute Z matrix: Z = W %*% J 
  Z_matrix <- compute_Z_matrix_generic(W, J, predictors, basis_functions)


  # Return the computed Y values
  list(W = W, Z = Z_matrix, Y = Y, J = J, B = Beta_matrix, U = U, X = X, basis_objs = basis_objs, basis_values = basis_values,beta_point_values = beta_point_values)
}


