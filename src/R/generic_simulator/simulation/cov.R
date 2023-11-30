library(MASS)
library(here)
library(fda)

source(here("src", "R",  "generic_simulator",   "utils" , "covariance_utilities.R"))

simulate_cov_data <- function(mu_funcs, cov_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr=c(FALSE,100)) {
  print(noise_snr)
  # Simulate functional features
  U <- simulate_functional_features(mu_funcs, cov_funcs, observations, time_domains)
  Betas <- create_beta_curves(beta_funcs, time_domains)
  
  X <- U
  # Add noise to Y (if specified)
  
  # X <- X + add_snr_noise( X, noise_snr[1])

  # Compute Y values
  Y <- compute_Y_values_with_func(U, Betas, observations, predictors, time_domains, intercept)$Y

  # Add noise to Y (if specified)
  Y <- Y + add_snr_noise(Y, noise_snr[2])


  return(list(U = NULL, X = X, Y = Y)) # U is not used in this simulation type
}





#' Simulate data
#' This function generates simulated data based on a mean function and a covariance function.
#' 
#' @param mu_func a function that takes in a time domain and returns the mean function evaluated at each time point.
#' @param cov_func_data a list containing the covariance function parameters: sig2, rho, and decay_type.
#' @param time_domain a vector of time points at which to evaluate the mean and covariance functions.
#' @param n the number of samples to generate.
#' 
#' @return A matrix of simulated data with n rows and length(time_domain) columns.
simulate_data <- function(mu_func, cov_func_data, time_domain, n) {
  mu <- mu_func(time_domain)
  # extract sig2, rho, decay_type from cov_func_data
  sig2 <- cov_func_data[[1]]
  rho <- cov_func_data[[2]]
  decay_type <- cov_func_data[[3]]
  
  Sigma <- generate_covariance_function(sig2, rho, decay_type)(time_domain, time_domain)
  MASS::mvrnorm(n, mu, Sigma)
}


#' Simulate functional features
#'
#' This function simulates functional features based on the provided mean and covariance functions.
#' 
#' @param mu_funcs A list of mean functions.
#' @param cov_funcs_data A list of covariance functions.
#' @param n The number of observations to simulate.
#' @param time_domains A list of time domains for each functional feature.
#' @return A 3D array of simulated functional features.
#' 
#' @examples
#' # Simulate 100 observations of 2 functional features with mean functions mu1 and mu2, and covariance functions cov1 and cov2
#' time_domains <- list(seq(0, 1, length.out = 10), seq(0, 1, length.out = 20))
#' mu_funcs <- list(mu1, mu2)
#' cov_funcs_data <- list(cov1, cov2)
#' simulate_functional_features(mu_funcs, cov_funcs_data, 100, time_domains)
#' 
simulate_functional_features <- function(mu_funcs, cov_funcs_data, observations, time_domains) {
  X_list = Map(simulate_data, mu_funcs, cov_funcs_data, time_domains, MoreArgs = list(n = observations))

  predictors = length(mu_funcs)
  X <- array(0, dim = c(observations, predictors,  length(time_domains[[1]])))

  # Fill the 3D array with data from simulated_features
  for (i in 1:predictors) {
    X[, i, ] <- X_list[[i]]
  }
  return(X)
}