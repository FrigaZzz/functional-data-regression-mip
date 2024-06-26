library(MASS)
library(here)
library(fda)


simulate_paper_data <- function(mu_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr=c(FALSE, FALSE)){
  # Simulate true predictors
  U <- simulate_true_predictors_Ut(mu_funcs, observations, time_domains)
  Betas <- create_beta_curves(beta_funcs, time_domains)

  X <- U
  # Assuming `data` is your dataset
  X <- simulate_observations_Xt(X)



  # Compute Y values
  Y_init <- compute_Y_values_with_func(U, Betas, observations, predictors, time_domains, intercept)$Y
  Y = Y_init
  # Add noise to Y (if specified)
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
  


# Loop over each predictor to calculate the range for all observations of that predictor
for (m in 1:predictors) {
  # Extract the functional values for predictor m across all observations
  u_values_all_obs <- X[, m, ]
  
  # Calculate the range for this specific predictor m across all observations
  RY <- max(u_values_all_obs) - min(u_values_all_obs)
  
  # Loop over each observation for this predictor
  for (i in 1:observations) {
    # Retrieve the functional values for predictor m and observation i
    u_values <- X[i, m, ]
    
    # Calculate the standard deviation of the error term based on the range
    sd_epsilon <-(0.025* RY)
    
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
compute_amplitude_norm <- function(Y, error_sd = 0.05) {
  observations <- length(Y)
  
  # Calculate the range of Y
  Y_range <- max(Y) - min(Y)
  
  # Calculate the standard deviation of the error term based on the range
  sd_epsilon <- error_sd * Y_range
  
  # Generate the error term vector
  error_vector <- rnorm(observations, mean = 0, sd = sd_epsilon)
  
  return(error_vector)
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

