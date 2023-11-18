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
simulate_true_predictors_Ut <- function(mu_funcs, n, time_domains,coef_list) {
  predictors = length(mu_funcs)
  measurements = length(time_domains[[1]])
  X <- array(0, dim = c(n, predictors, measurements))
  # Generate the functional data for each observation and each predictor
  for (i in 1:n) {
    for (j in 1:predictors) {
      # Extract the coefficients for the j-th predictor of the i-th observation
      current_coefs <- lapply(coef_list[[j]], function(row) row[i])
      # Apply the mu_funcs function with the time domain for the j-th predictor
      # and the current coefficients for the i-th observation
      X[i, j, ] <- mu_funcs[[j]](time_domains[[j]], current_coefs)
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
simulate_observations_Xt <- function(X) {
  measurements <- dim(X)[3]
  predictors <- dim(X)[2]
  observations <- dim(X)[1]
  # Initialize the array for the normalized functional covariates
  FX <- array(0, dim = c(observations, predictors, measurements))

  for (i in 1:observations) {
    for (m in 1:predictors) {
      # Retrieve the functional values for predictor m and observation i
      u_values <- X[i, m, ]
      
      # Calculate the range for this set of functional values
      RY <- max(u_values) - min(u_values)
      
      # Calculate the standard deviation of the error term based on the range
      sd_epsilon <- (0.025 * RY)
      
      # Generate the error term for this observation and predictor
      epsilon_im <- rnorm(measurements, mean = 0, sd = sd_epsilon)
      
      # Add the error term to the functional covariate
      FX[i, m, ] <- u_values + epsilon_im
    }
  }
  
  return(FX)
}


#' Compute Y matrix given observations, predictors, basis functions, B matrix, beta_0, and Z matrix
#'
#' @param observations Number of observations
#' @param predictors Number of predictors
#' @param basis_functions Basis functions
#' @param B_matrix B matrix
#' @param beta_0 Intercept
#' @param Z_matrix Z matrix
#'
#' @return Y matrix
#'
#' @examples
#' compute_Y_values(10, 5, basis_functions, B_matrix, 0, Z_matrix)
compute_Y_values <- function(mu_funcs, beta_funcs, observations, predictors, time_domains, intercept,coef_list=NULL) {
  Y <- numeric(observations)
  
  for (i in 1:observations) {
    for (j in 1:predictors) {
      # Define a function that is the product of the functional predictor and the coefficient function
      current_coefs <- lapply(coef_list[[j]], function(row) row[i])

      integrand <- function(t) {
        mu = mu_funcs[[j]](t,current_coefs)
        beta = beta_funcs[[j]](t)
        mu * beta
        # Perform numerical integration of the integrand over the time domain Tm
      
        }
        Y[i] <- Y[i] + integrate(integrand, lower = min(time_domains[[j]]),   upper = max(time_domains[[j]]))$value
    }
  }
  
  return(list(Y=Y))
}



#' Compute amplitude normalization factor for functional data regression
#'
#' This function computes the amplitude normalization factor for functional data regression.
#' The function takes in a vector of observations Y and an error standard deviation error_sd.
#' It returns a vector of random errors Epsilon with the same length as Y.
#'
#' @param Y A vector of observations
#' @param error_sd The standard deviation of the error term
#' 
#' @return A vector of random errors Epsilon with the same length as Y
#'
#' @examples
#' Y <- rnorm(100)
#' compute_amplitude_norm(Y, error_sd = 0.05)
#'
compute_amplitude_norm<- function(Y, error_sd = 0.05) {
  observations <- length(Y)
  Epsilon <- numeric(observations)
  Ry <- max(Y) - min(Y)  # Calculate the range of Y
  Epsilon <- rnorm(observations, mean = 0, sd = (error_sd * Ry))
  return(Epsilon)
}

#' Compute noise for a given set of observations
#'
#' This function generates random noise for a given set of observations.
#' The noise is generated from a normal distribution with mean 0 and standard deviation error_sd.
#'
#' @param Y A numeric vector of observations
#' @param error_sd A numeric value representing the standard deviation of the error term
#'
#' @return A numeric vector of random noise
#'
#' @examples
#' Y <- c(1, 2, 3, 4, 5)
#' compute_noise(Y, error_sd = 0.1)
#'
compute_noise <- function(Y, error_sd = 0.05) {
  observations <- length(Y)
  noise <- rnorm(observations, mean = 0, sd = error_sd)
  return(noise)
}