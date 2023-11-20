library(MASS)
library(here)
library(fda)

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


#' 
#' This function simulates functional features for the paper using the provided mean functions,  number of observations and time domains.
#' 
#' @param mu_funcs A list of mean functions.
#' @param observations The number of observations.
#' @param time_domains A list of time domains.
#' 
#' @return An array of simulated functional features.
#' 
simulate_true_predictors_Ut <- function(mu_funcs, observations, time_domains) {
  predictors = length(mu_funcs)
  measurements = length(time_domains[[1]])
  coef_list <- list(
    '1' = list(a1 = rnorm(observations, mean = -4, sd = 3), a2 = rnorm(observations, mean = 7, sd = 1.5)),
    '2' = list(b1 = runif(observations, min = 3, max = 7), b2 = rnorm(observations, mean = 0, sd = 1)),
    '3' = list(c1 = rnorm(observations, mean = -3, sd = sqrt(1.2^2)), c2 = rnorm(observations, mean = 2, sd = sqrt(0.5^2)), c3 = rnorm(observations, mean = -2, sd = 1)),
    '4' = list(d1 = rnorm(observations, mean = -2, sd = 1), d2 = rnorm(observations, mean = 3, sd = sqrt(1.5^2))),
    '5' = list(e1 = runif(observations, min = 2, max = 7), e2 = rnorm(observations, mean = 2, sd = sqrt(0.4^2))),
    '6' = list(f1 = rnorm(observations, mean = 4, sd = sqrt(2^2)), f2 = rnorm(observations, mean = -3, sd = sqrt(0.5^2)), f3 = rnorm(observations, mean = 1, sd = 1))
  )

  X <- array(0, dim = c(observations, predictors, measurements))
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


#' Create Beta Curves
#'
#' Generates beta curves for functional data based on provided functions and time domains.
#'
#' @param beta_funcs Functions to generate beta values.
#' @param time_domains Time domains for beta functions.
#' @return Array of beta curves.
create_beta_curves <- function(beta_funcs, time_domains) {
  predictors = length(beta_funcs)
  Betas <- array(0, dim = c(predictors, length(time_domains[[1]])))
  for (j in 1:predictors) {
    Betas[j, ] <- beta_funcs[[j]](time_domains[[j]])
  }
  return(Betas)
}


#' Apply amplitude normalization to predictor data and add random error terms
#'
#' This function applies amplitude normalization to the predictor data and adds random error terms to it.
#' 
#' @param X array of predictor data
#' @param observations number of observations
#' @param predictors number of predictors
#' @return An array of predictor data with added error terms.
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

#' Add SNR Noise to Signal
#'
#' Adds noise to a signal to achieve a specified Signal-to-Noise Ratio (SNR).
#'
#' @param Y Signal to which noise will be added.
#' @param snr_db Desired SNR in decibels.
#' @return Signal with added noise.
apply_snr_to_X <- function(X, snr_db = 100) {
  # Determine the size of X
  num_observations <- dim(X)[1]
  num_predictors <- dim(X)[2]
  num_time_points <- dim(X)[3]

  # Iterate over each predictor
  for (predictor in 1:num_predictors) {
    for (time_point in 1:num_time_points) {
      # Extract the measurements for this predictor across all observations at the given time point
      measurements <- X[, predictor, time_point]

      # Apply the SNR function to these measurements
      measurements_noisy <- add_snr_noise(measurements, snr_db)

      # Update the measurements in X
      X[, predictor, time_point] <- X[, predictor, time_point] + measurements_noisy
    }
  }

  return(X)
}

#' Compute Y Values
#'
#' Computes Y values (response variable) using functional data, beta curves, and other parameters.
#'
#' @param X_data Functional data for predictors.
#' @param beta_curves Beta curves for predictors.
#' @param observations Number of observations.
#' @param predictors Number of predictors.
#' @param time_domains Time domains for predictors.
#' @param intercept Intercept value for Y computation.
#' @return List containing computed Y values.
compute_Y_values <- function(X_data, beta_curves, observations, predictors, time_domains, intercept = 0) {
  Y <- numeric(observations)
  basis_functions <- 6

  # Create basis functions for each predictor
  basis_list <- lapply(time_domains, function(td) {
    create.bspline.basis(rangeval = c(min(td), max(td)), nbasis = basis_functions)
  })

  # Smooth beta_curves
  beta_fd_list <- lapply(1:predictors, function(j) {
    smooth.basis(time_domains[[j]], beta_curves[j, ], basis_list[[j]])$fd
  })

  # Compute Y values for each observation
  for (i in 1:observations) {
    Y[i] <- intercept
    for (j in 1:predictors) {
      # Extract the time series for the current observation and predictor
      current_series <- X_data[i, j, ]

      # Smooth the current series into a functional data object
      X_fd <- smooth.basis(time_domains[[j]], current_series, basis_list[[j]])$fd

      # Compute the product of X_fd and beta_fd as a new functional data object
      product_fd <- X_fd * beta_fd_list[[j]]

      # Create a constant function for integration
      constant_fd <- fd(matrix(1, nrow = basis_functions, ncol = 1), basis_list[[j]])

      # Perform functional integration over the time domain
      Y[i] <- Y[i] + inprod(product_fd, constant_fd)
    }
  }

  return(list(Y = Y))
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
#'
compute_noise <- function(Y, error_sd = 0.05) {
  observations <- length(Y)
  noise <- rnorm(observations, mean = 0, sd = error_sd)
  return(noise)
}




#' Compute Noise and Add to Signal
#'
#' This function adds noise to a given signal `Y` such that the resulting noisy signal has an approximate Signal-to-Noise Ratio (SNR) specified by `snr_db`. 
#' The noise added is normally distributed.
#'
#' @param Y Numeric vector representing the original signal to which noise will be added.
#' @param snr_db Desired Signal-to-Noise Ratio in decibels (dB). This parameter specifies the relative power of the signal to the noise.
#' @param error_sd Optional parameter, standard deviation of the error, with a default value of 0.05. This parameter is used to adjust the scale of the noise if needed.
#'
#' @return Returns the original signal `Y` with added noise. The resulting noisy signal has an SNR approximately equal to the specified `snr_db`.
#'
add_snr_noise <- function(Y, snr_db = 0) {
  # Calculate the power of the signal
  signal_power <- var(Y)

  # Convert SNR from dB to linear scale
  snr_linear <- 10^(snr_db / 10)

  # Calculate the required noise power
  noise_power <- signal_power / snr_linear

  # Calculate the standard deviation for the noise
  noise_sd <- sqrt(noise_power)
  # Generate noise
  Epsilon <- rnorm(length(Y), mean = 0, sd = noise_sd)


  return(Epsilon)
}
