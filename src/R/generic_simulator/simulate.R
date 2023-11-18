
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
#' @param error_sd A numeric value for the standard deviation of the error terms
#' @param seed An integer value for the random seed
#' @param noise_sd A numeric value for the standard deviation of the amplitude normalization
#'
#' @return A list containing the computed Y values

generate_data <- function(
    predictors, observations, measurements , basis_functions, intercept = 0, norder,
    mu_funcs, beta_funcs, time_domains ,cov_funcs ,error_sd =0, seed = 2000, noise_sd = 0.05) {

  # Set seed for reproducibility
  set.seed(seed)

  coef_list <- list(
      '1' = list(a1 = rnorm(observations, mean = -4, sd = 3), a2 = rnorm(observations, mean = 7, sd = 1.5)),
      '2' = list(b1 = runif(observations, min = 3, max = 7), b2 = rnorm(observations, mean = 0, sd = 1)),
      '3' = list(c1 = rnorm(observations, mean = -3, sd = sqrt(1.2^2)), c2 = rnorm(observations, mean = 2, sd = sqrt(0.5^2)), c3 = rnorm(observations, mean = -2, sd = 1)),
      '4' = list(d1 = rnorm(observations, mean = -2, sd = 1), d2 = rnorm(observations, mean = 3, sd = sqrt(1.5^2))),
      '5' = list(e1 = runif(observations, min = 2, max = 7), e2 = rnorm(observations, mean = 2, sd = sqrt(0.4^2))),
      '6' = list(f1 = rnorm(observations, mean = 4, sd = sqrt(2^2)), f2 = rnorm(observations, mean = -3, sd = sqrt(0.5^2)), f3 = rnorm(observations, mean = 1, sd = 1))
  )
  # 1. Simulate functional features
  U <- NULL
  X <- NULL
  # Pre-generate all coefficients for each predictor and each observation

  U <- simulate_true_predictors_Ut(mu_funcs, observations, time_domains,coef_list) 
  # Apply amplitude normalization
  X <- simulate_observations_Xt(U)

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
  W <- compute_W_matrix_generic(X, basis_functions, time_domains, basis_objs)

  # 6. Compute Z matrix: Z = W %*% J 
  Z_matrix <- compute_Z_matrix_generic(W, J, predictors, basis_functions)

  # 7. Compute Y values
  Y <- compute_Y_values(mu_funcs, beta_funcs, observations, predictors, time_domains,intercept,coef_list)$Y
  # Gu <- result_Y$Gu

  # 8. Apply error terms to Y values

  Y <- Y + compute_amplitude_norm(Y)

  # Return the computed Y values
  list(W = W, Z = Z_matrix, Y = Y, J = J, B = Beta_matrix, U = U, X = X, basis_objs = basis_objs, basis_values = basis_values,beta_point_values = beta_point_values)
}


