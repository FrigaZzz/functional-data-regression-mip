library(MASS)
library(here)
library(fda)




simulate_paper_data <- function(mu_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr=c(FALSE, FALSE)){
  # Simulate true predictors
  U <- simulate_true_predictors_Ut(mu_funcs, observations, time_domains)
  Betas <- create_beta_curves(beta_funcs, time_domains)

  X <- U
  if(noise_snr[1] == TRUE)
    X <- simulate_observations_Xt(X)

  # Compute Y values
  Y <- compute_Y_values_generic(U, Betas, observations, predictors, time_domains, intercept)$Y

  # Add noise to Y (if specified)
  if(noise_snr[1] == TRUE){
    Y <- Y + compute_amplitude_norm(Y)
  }

  return(list(U = U, X = X, Y = Y))
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
  for (i in 1:observations) {
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
      sd_epsilon <- (0.005 * RY)
      # Generate the error term for this observation and predictor
      epsilon_im <- rnorm(measurements, mean = 0, sd = sd_epsilon)
      
      # Add the error term to the functional covariate
      FX[i, m, ] <- u_values + epsilon_im
    }
  }
  
  return(FX)
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
#' compute_amplitude_norm(Y, error_sd = 1)
#'
compute_amplitude_norm<- function(Y, error_sd = 1) {
  observations <- length(Y)
  Epsilon <- numeric(observations)
  Ry <- max(Y) - min(Y)  # Calculate the range of Y
  Epsilon <- rnorm(observations, mean = 0, sd = (error_sd * Ry))
  return(Epsilon)
}