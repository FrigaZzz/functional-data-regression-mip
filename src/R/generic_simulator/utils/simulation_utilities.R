library(MASS)
library(here)

source(here("src", "R",  "generic_simulator",   "utils" , "covariance_utilities.R"))

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
#' 
#' @return A 3D array of simulated functional features.
#' 
#' @examples
#' # Simulate 100 observations of 2 functional features with mean functions mu1 and mu2, and covariance functions cov1 and cov2
#' time_domains <- list(seq(0, 1, length.out = 10), seq(0, 1, length.out = 20))
#' mu_funcs <- list(mu1, mu2)
#' cov_funcs_data <- list(cov1, cov2)
#' simulate_functional_features(mu_funcs, cov_funcs_data, 100, time_domains)
#' 
simulate_functional_features <- function(mu_funcs, cov_funcs_data, n, time_domains) {
  X_list = Map(simulate_data, mu_funcs, cov_funcs_data, time_domains, MoreArgs = list(n = n))

  predictors = length(mu_funcs)
  X <- array(0, dim = c(n, predictors,  length(time_domains[[1]])))

  # Fill the 3D array with data from simulated_features
  for (i in 1:predictors) {
    X[, i, ] <- X_list[[i]]
  }
  return(X)
}


#' 
#' This function simulates functional features for the paper using the provided mean functions,  number of observations and time domains.
#' 
#' @param mu_funcs A list of mean functions.
#' @param n The number of observations.
#' @param time_domains A list of time domains.
#' 
#' @return An array of simulated functional features.
#' 
simulate_functional_features_paper <- function(mu_funcs, n, time_domains) {
  predictors = length(mu_funcs)
  X <- array(0, dim = c(n, predictors,  length(time_domains[[1]])))
  coef_list <- list(
      list(a1 = rnorm(1, mean = -4, sd = 3), a2 = rnorm(1, mean = 7, sd = 1.5)),
      list(b1 = runif(1, min = 3, max = 7), b2 = rnorm(1, mean = 0, sd = 1)),
      list(c1 = rnorm(1, mean = -3, sd = sqrt(1.2^2)), c2 = rnorm(1, mean = 2, sd = sqrt(0.5^2)), c3 = rnorm(1, mean =    -2, sd = 1)),
      list(d1 = rnorm(1, mean = -2, sd = 1), d2 = rnorm(1, mean = 3, sd = sqrt(1.5^2))),
      list(e1 = runif(1, min = 2, max = 7), e2 = rnorm(1, mean = 2, sd = sqrt(0.4^2))),
      list(f1 = rnorm(1, mean = 4, sd = sqrt(2^2)), f2 = rnorm(1, mean = -3, sd = sqrt(0.5^2)), f3 = rnorm(1, mean = 1,     sd = 1))
  )
  for (i in 1:n) {
    for (j in 1:predictors) {
      X[i, j, ] <- mu_funcs[[j]](time_domains[[j]])
    }
  }
  return(X)
}


#' Apply amplitude normalization to predictor data and add random error terms
#'
#' This function applies amplitude normalization to the predictor data and adds random error terms to it.
#' 
#' @param X array of predictor data
#' @param observations number of observations
#' @param predictors number of predictors
#' 
#' @return An array of predictor data with added error terms.
#'
#' @examples
#' apply_amplitude_norm(X, observations, predictors)
#'
apply_amplitude_norm <- function(X, observations, predictors) {
  measurements <- dim(X)[3]
  # Initialize the array for the normalized functional covariates
  FX <- array(NA, dim = c(observations, predictors, measurements))

  for (i in 1:observations) {
    for (m in 1:predictors) {
      # Retrieve the functional values for predictor m and observation i
      u_values <- X[i, m, ]
      
      # Calculate the range for this set of functional values
      RY <- max(u_values) - min(u_values)
      
      # Calculate the standard deviation of the error term based on the range
      sd_epsilon <- sqrt(0.025 * RY^2)
      
      # Generate the error term for this observation and predictor
      epsilon_im <- rnorm(measurements, mean = 0, sd = sd_epsilon)
      
      # Add the error term to the functional covariate
      FX[i, m, ] <- u_values + epsilon_im
    }
  }
  
  return(FX)
}




