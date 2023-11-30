library(MASS)
library(here)
library(fda)

source(here("src", "R",  "generic_simulator",   "utils" , "covariance_utilities.R"))
source(here("src", "R",  "generic_simulator",   "simulation" , "cov.R"))
source(here("src", "R",  "generic_simulator",   "simulation" , "paper.R"))




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
compute_Y_values_smoothed <- function(X_data, beta_curves, observations, predictors, time_domains, intercept = 0, coef_list = null) {
  Y <- numeric(observations)
  basis_functions <- 10

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
      current_predictor <- X_data[i, j, ]

      # Smooth the current series into a functional data object
      X_fd <- smooth.basis(time_domains[[j]], current_predictor, basis_list[[j]])$fd

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



compute_Y_values_with_func <- function(X_data, beta_curves, observations, predictors, time_domains, intercept = 0) {
  Y <- numeric(observations)

  for (i in 1:observations) {
    Y[i] <- intercept
    for (j in 1:predictors) {
      current_series <- X_data[i, j, ]
      beta_series <- beta_curves[j, ]
      td <- time_domains[[j]]
      prod = current_series * beta_series
      area = trapz(td, prod)
      Y[i] <- Y[i] + area
    }
  }

  return(list(Y = Y))
}

# when times in time domains are not equidistant
compute_Y_values_generic_ <- function(X_data, beta_curves, observations, predictors, time_domains, intercept = 0) {
  # Initialize the response vector
  Y <- numeric(observations)
  
  # Loop through each observation
  for (i in 1:observations) {
    # Start with the intercept for each Y[i]
    Y[i] <- intercept
    
    # Loop through each predictor
    for (j in 1:predictors) {
      # Extract the time series and beta values for the current observation and predictor
      current_series <- X_data[i, j, ]
      beta_series <- beta_curves[j, ]
      
      # Extract the time domain for the current predictor
      td <- time_domains[[j]]
      
      # Initialize integral approximation for this predictor
      integral_approximation <- 0
      
      # Loop through the time domain points
      for (k in 1:(length(td) - 1)) {
        # Calculate the non-uniform width of the current trapezoid
        h_k <- td[k + 1] - td[k]
        
        # Calculate the average height of the current trapezoid
        average_height <- (current_series[k] * beta_series[k] + current_series[k + 1] * beta_series[k + 1]) / 2
        
        # Add the area of the current trapezoid to the integral approximation
        integral_approximation <- integral_approximation + h_k * average_height
      }
      
      # Update Y[i] with the contribution from the current predictor
      Y[i] <- Y[i] + integral_approximation
    }
  }
  
  # Return the computed Y values
  return(list(Y = Y))
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
add_snr_noise <- function(Y, snr_linear = 5) {
  # Calculate the power of the signal
  signal_power <- var(Y)

  # Calculate the required noise power
  noise_power <- signal_power / snr_linear

  # Calculate the standard deviation for the noise
  noise_sd <- sqrt(noise_power)
  # Generate noise
  Epsilon <- rnorm(length(Y), mean = 0, sd = noise_sd)


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



trapz <- function(x, y) {
    if (missing(y)) {
        if (length(x) == 0) return(0)
        y <- x
        x <- seq(along=x)
    }
    if (length(x) == 0 && length(y) == 0) return(0)
    if (!(is.numeric(x) || is.complex(x)) ||
            !(is.numeric(y) || is.complex(y)) )
        stop("Arguments 'x' and 'y' must be real or complex vectors.")
    m <- length(x)
    if (length(y) != m)
        stop("Arguments 'x', 'y' must be vectors of the same length.")
    if (m <= 1) return(0.0)

    # z <- sum((x[2:m] - x[1:(m-1)]) * (y[1:(m-1)] + y[2:m]))
    # return(0.5 * z)

    xp <- c(x, x[m:1])
    yp <- c(numeric(m), y[m:1])
    n <- 2*m
    p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1]
    p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n]

    return(0.5*(p1-p2))
}
