# Load necessary libraries
library(fda)

#' Create a basis object for a single functional predictor
#'
#' @param basis_functions The number of basis functions to use
#' @param bm_range The time domain for the functional predictor
#' @param norder The order of the B-spline basis functions
#'
#' @return A basis object for the functional predictor
create_1_basis <- function(basis_functions, bm_range, norder) {
  # Loop over each functional predictor
    # Use the specific time domain for each functional predictor to create  the basis object
    basis_obj<- fda::create.bspline.basis(rangeval = c(min (bm_range), max(bm_range)), nbasis = basis_functions, norder = norder)
  return(basis_obj)
}

#' Create a list of basis objects for multiple functional predictors
#'
#' @param basis_functions The number of basis functions to use
#' @param time_domains A list of time domains for each functional predictor
#' @param norder The order of the B-spline basis functions
#' @param predictors The number of functional predictors
#'
#' @return A list of basis objects for each functional predictor
create_basis <- function(basis_functions, time_domains, norder, predictors) {
  basis_objs <- vector("list", length = predictors)
  # Loop over each functional predictor
  for (m in 1:predictors) {
    # Use the specific time domain for each functional predictor to create  the basis object
    bm_range <- time_domains[[m]]
    basis_objs[[m]] <- fda::create.bspline.basis(rangeval = c(min (bm_range), max(bm_range)), nbasis = basis_functions, norder = norder)
  }
  return(basis_objs)
}

#' Function to expand functional data into B-spline basis
#' 
#' @param data_matrix A matrix of functional data
#' @param time A vector of time points
#' @param basis_obj A B-spline basis object
#' 
#' @return A matrix of coefficients for the B-spline basis
expand_functional_data <- function(data_matrix, time, basis_obj) {
  fd_list <- apply(data_matrix, 1, function(row) fda::smooth.basis(time, row, basis_obj)$fd)
  do.call(rbind, lapply(fd_list, function(fd) t(fd$coefs)))
}

#' Smooth beta coefficients into the B-spline basis
#' 
#' @param beta_func A function that returns beta values at given time points
#' @param time A vector of time points
#' @param basis_obj A B-spline basis object
#' 
#' @return A vector of smoothed beta coefficients
smooth_beta_function <- function(beta_func,time,basis_obj) {
  beta_values <- beta_func(time)
  fdPar_obj <- fdPar(basis_obj)
  smoothed_beta <- smooth.basis(time, beta_values, fdPar_obj)
  beta_coefs <- coef(smoothed_beta$fd)
  beta_coefs <- beta_coefs[, 1]
  return(beta_coefs)
}



#' Smooth beta coefficients into the B-spline basis for multiple predictors with different time domains
#' 
#' @param beta_funcs A list of functions that return beta values at given time points
#' @param num_basis The number of basis functions
#' @param time_domains A list of vectors of time points for each predictor
#' @param basis_objs A list of B-spline basis objects for each predictor
#' 
#' @return A list containing the smoothed beta coefficients, the beta point values, and the basis values
smooth_betas_generic <- function(beta_funcs,num_basis, time_domains, basis_objs) {
  
  num_predictors <- length(beta_funcs)
  measurements <- length(time_domains[[1]])
  # Initialize a matrix to store the smoothed beta coefficients
  beta_matrix <-  array(0, dim = c(num_predictors,num_basis))
  beta_point_values <-  array(0, dim = c(num_predictors, measurements))
  basis_values <-  array(0, dim = c(num_predictors,measurements,num_basis))
  for (i in 1:num_predictors) {
    beta_func <- beta_funcs[[i]]
    basis_obj <- basis_objs[[i]]
    time <- time_domains[[i]]

    beta_values <- beta_func(time)
    if(length(beta_values) == 1){ # if beta is a 0
      beta_values <- rep(0, length(time))
    }
    beta_point_values[i,] <- beta_values
    basis_values[i,,] <- eval.basis(time, basis_obj)
    fdPar_obj <- fda::fdPar(basis_obj)
    smoothed_beta <- fda::smooth.basis(time, beta_values, fdPar_obj)
    beta_matrix[i,] <- smoothed_beta$fd$coefs
  }
  return(list(beta_matrix = beta_matrix, beta_point_values = beta_point_values, basis_values = basis_values))

}


test = FALSE
if (test) {
  # Generate some example data
  n_obs <- 100
  n_measurements <- 100
  n_funcs <- 5
  time <- seq(0, 1, length.out = n_measurements)

  # Create data matrix with dimensions consistent with time vector
  data_matrix <- matrix(rnorm(n_obs * n_measurements), nrow = n_obs, ncol = n_measurements)

  # Define B-spline basis
  n_basis <- 10
  bspline_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis, norder = 4)

  # Expand data into B-spline basis
  expanded_data <- expand_functional_data(data_matrix, time, bspline_basis)

  # Define beta function
  beta_func <- function(t) sin(2 * pi * t)

  # Smooth beta function into B-spline basis
  beta_coefs <- smooth_beta_function(beta_func, time, bspline_basis)

}
