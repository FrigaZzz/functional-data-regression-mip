library(MASS)
library(here)
library(fda)


simulate_paper_data <- function(mu_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr=c(FALSE, FALSE)){
  # Simulate true predictors
  U <- simulate_true_predictors_Ut(mu_funcs, observations, time_domains)
  Betas <- create_beta_curves(beta_funcs, time_domains)

  X <- U

  X <- simulate_observations_Xt(X)

  # Compute Y values
  Y_init <- compute_Y_values_with_func(U, Betas, observations, predictors, time_domains, intercept)$Y


  Y <- Y_init + compute_amplitude_norm(Y_init)
  print(compute_snr(U, X))
  print(compute_snr(Y_init, Y))
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
  coef_list <- generate_coefficients(observations, coef_specs)
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

  # for (i in 1:observations) {
  #   for (m in 1:predictors) {
  #     # Retrieve the functional values for predictor m and observation i
  #     u_values <- X[i, m, ]
      
  #     # Calculate the range for this set of functional values
  #     RY <- max(u_values) - min(u_values)
      
  #     # Calculate the standard deviation of the error term based on the range
  #     sd_epsilon <- (0.025 * RY)
  #     # Generate the error term for this observation and predictor
  #     epsilon_im <- rnorm(measurements, mean = 0, sd = sd_epsilon)
      
  #     # Add the error term to the functional covariate
  #     FX[i, m, ] <- u_values + epsilon_im
  #   }
  # }
  for (m in 1:predictors) {
    # Initialize variables to track the min and max values
    min_value = Inf
    max_value = -Inf
    
    # Loop through all observations to find the min and max values of u_alpha_m(t)
    for (i in 1:observations) {
        # Retrieve the functional values for predictor m and observation i
        u_values <- X[i, m, ]
        
        # Update min and max values if necessary
        min_value = min(min(u_values), min_value)
        max_value = max(max(u_values), max_value)
    }
    
    # Calculate the range for this predictor across all observations
    RY = max_value - min_value
    
    # Calculate the standard deviation of the error term based on the range
    sd_epsilon = 0.025 * RY
    
    # Generate the error term for this predictor
    
    # Apply the error term to all observations of F(X) for this predictor
    for (i in 1:observations) {
        # Retrieve the functional values for predictor m and observation i
        u_values = X[i, m, ]
          
        epsilon_im <- rnorm(measurements, mean = 0, sd = sd_epsilon)

        # Add the error term to the functional covariate
        FX[i, m,] = u_values + epsilon_im
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
compute_amplitude_norm<- function(Y, error_sd = 0.05) {
  observations <- length(Y)
  Epsilon <- numeric(observations)
  Ry <- max(Y) - min(Y)  # Calculate the range of Y
  Epsilon <- rnorm(observations, mean = 0, sd = (error_sd * Ry))
  return(Epsilon)
}


generate_coefficients <- function(observations, coef_specs) {
  coef_list <- list()
  for (predictor in names(coef_specs)) {
    predictor_specs <- coef_specs[[predictor]]
    predictor_coefs <- list()
    for (coef_name in names(predictor_specs)) {
      spec <- predictor_specs[[coef_name]]
      if (spec$type == "norm") {
        predictor_coefs[[coef_name]] <- rnorm(observations, mean = spec$mean, sd = spec$sd)
      } else if (spec$type == "unif") {
        predictor_coefs[[coef_name]] <- runif(observations, min = spec$min, max = spec$max)
      }
    }
    coef_list[[predictor]] <- predictor_coefs
  }
  return(coef_list)
}

