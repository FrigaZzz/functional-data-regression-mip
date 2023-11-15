library(fda)

# Compute W matrix
# This function computes the B-spline coefficients for each predictor in the input matrix FX
# and returns a 3D array containing the coefficients for each observation and predictor.
# Args:
#   FX: input matrix of predictors
#   num_basis: number of B-spline basis functions
#   time: vector of time points
#   basis_obj: B-spline basis object
# Returns:
#   A 3D array containing the B-spline coefficients for each observation and predictor.
#' @param FX input matrix of predictors
#' @param num_basis number of B-spline basis functions
#' @param time vector of time points
#' @param basis_obj B-spline basis object
#' @return A 3D array containing the B-spline coefficients for each observation and predictor.
compute_W_matrix <- function(FX,num_basis, time, basis_obj) {
  num_predictors <- dim(FX)[2]
  num_observations <- dim(FX)[1]
  num_basis <- num_basis
  
  # Initialize a 3D array to store the B-spline coefficients
  expanded_data_array <-  array(0, dim = c(num_observations,num_predictors,num_basis))

  for (i in 1:num_observations) {
    for (j in 1:num_predictors) {
      data_vector <- FX[i, j, ]
      fd_obj <- fda::smooth.basis(time, data_vector, basis_obj)$fd
      expanded_data_array[i,j,] <- fd_obj$coefs
    }
  }
  
  return(expanded_data_array)
}

#' Compute W matrix, Different time domains
#'
#' This function computes the W matrix for functional data regression using B-splines. 
#' The W matrix is a 3D array of B-spline coefficients for each observation and predictor.
#' 
#' @param FX A 3D array of functional data predictors
#' @param num_basis The number of B-spline basis functions to use
#' @param time_domains A list of time domains for each predictor
#' @param basis_obj_list A list of B-spline basis objects for each predictor
#' 
#' @return A 3D array of B-spline coefficients for each observation and predictor
#'
#' @examples
#' # Generate sample data
#' time_domains <- list(seq(0, 1, length.out = 10), seq(0, 1, length.out = 20))
#' basis_obj_list <- lapply(time_domains, fda::create.bspline.basis, degree = 3, norder = 4)
#' FX <- array(rnorm(200), dim = c(10, 2, 10))
#' 
#' # Compute W matrix
#' W <- compute_W_matrix_paper(FX, num_basis = 10, time_domains, basis_obj_list)
#' 
compute_W_matrix_paper <- function(FX, num_basis, time_domains, basis_obj_list) {
  num_predictors <- dim(FX)[2]
  num_observations <- dim(FX)[1]
  num_basis <- num_basis
  
  # Initialize a 3D array to store the B-spline coefficients
  expanded_data_array <-  array(0, dim = c(num_observations,num_predictors,num_basis))

  for (i in 1:num_observations) {
    for (j in 1:num_predictors) {
      data_vector <- FX[i, j, ]
      time <- time_domains[[j]]
      basis_obj <- basis_obj_list[[j]]      
      fd_obj <- fda::smooth.basis(time, data_vector, basis_obj)$fd
      expanded_data_array[i,j,] <- fd_obj$coefs
    }
  }
  
  return(expanded_data_array)
}


#' Compute the J matrix (inner products for basis functions)
#'
#' This function computes the J matrix, which represents the inner products for basis functions.
#' 
#' @param basis_obj an object of class \code{fd} containing the basis functions
#' 
#' @return a matrix representing the inner products for basis functions
#' 
#' @examples
#' # Load the required packages
#' library(fda)
#' 
#' # Create a basis object
#' basis_obj <- create.bspline.basis(rangeval = c(0, 1), nbasis = 10, norder = 4)
#' 
#' # Compute the J matrix
#' J_matrix <- compute_J_matrix(basis_obj)
#' 
# Compute the J matrix (inner products for basis functions)
compute_J_matrix <- function(basis_obj) {
  fda::inprod(basis_obj, basis_obj)
}


#' Compute Z matrix
#'
#' This function computes the Z matrix for functional data regression models.
#' The Z matrix is a three-dimensional array with dimensions observations x predictors x basis_functions.
#' Each element of the Z matrix is the product of the corresponding element of the W_array and J.
#'
#' @param W_array A three-dimensional array with dimensions observations x predictors x basis_functions.
#' @param J A matrix with dimensions basis_functions x number_of_functions.
#' @param predictors An integer indicating the number of predictors.
#' @param basis_functions An integer indicating the number of basis functions.
#'
#' @return A three-dimensional array with dimensions observations x predictors x basis_functions.
#'
#' @examples
#' W_array <- array(1:12, dim = c(2, 2, 3))
#' J <- matrix(1:6, nrow = 3)
#' compute_Z_matrix(W_array, J, 2, 3)
#'
compute_Z_matrix <- function(W_array, J, predictors, basis_functions) {
  observations <- dim(W_array)[1]
  
  # Initializing Z_matrix 
  Z_matrix <- matrix(0, nrow = observations, ncol = predictors * basis_functions)

  # Looping through each observation
  for (i in 1:observations) {
    Z_matrix[i, ] <- W_array[i, , ] %*% J[ , ] 
  }
  Z_matrix <- array(Z_matrix, dim = c(observations,predictors,basis_functions))

  return(Z_matrix)
}





# Compute Y matrix
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
compute_Y_values <- function(observations, predictors, basis_functions, B_matrix, beta_0, Z_matrix) {
  Y <- rep(beta_0, observations)

  # Compute the sum of products of Z_matrix elements and B_matrix coefficients
  for (i in 1:observations) {
    Y[i] <-  sum(Z_matrix[i,,] * B_matrix)
  }

  return(Y)
}

#' Compute W matrix, Different time domains
#'
#' This function computes the W matrix for a given set of functional predictors and basis functions.
#' The W matrix is a 3D array that stores the B-spline coefficients for each predictor and observation.
#'
#' @param FX A 3D array of functional predictors
#' @param num_basis The number of basis functions to use
#' @param time_domains A list of time domains for each predictor
#' @param basis_obj_list A list of basis function objects for each predictor
#'
#' @return A 3D array of B-spline coefficients for each predictor and observation
#'
#' @examples
#' # Generate some sample data
#' time_domains <- list(seq(0, 1, length.out = 10), seq(0, 1, length.out = 20))
#' basis_obj_list <- lapply(time_domains, fda::create.bspline.basis, degree = 3, nbreaks = 6)
#' FX <- array(rnorm(200), dim = c(10, 2, 10))
#' 
#' # Compute the W matrix
#' W <- compute_W_matrix_generic(FX, num_basis = 6, time_domains, basis_obj_list)
#'
compute_W_matrix_generic <- function(FX, num_basis, time_domains, basis_obj_list) {
  num_predictors <- dim(FX)[2]
  num_observations <- dim(FX)[1]
  num_basis <- num_basis
  
  # Initialize a 3D array to store the B-spline coefficients
  expanded_data_array <-  array(0, dim = c(num_observations,num_predictors,num_basis))

  for (i in 1:num_observations) {
    for (j in 1:num_predictors) {
      data_vector <- FX[i, j, ]
      time <- time_domains[[j]]
      basis_obj <- basis_obj_list[[j]]      
      fd_obj <- fda::smooth.basis(time, data_vector, basis_obj)$fd
      expanded_data_array[i,j,] <- fd_obj$coefs
    }
  }
  
  return(expanded_data_array)
}



#' Compute the J matrix for functional data regression
#'
#' This function computes the J matrix for functional data regression, which is used in the estimation of the regression coefficients.
#' The J matrix is a three-dimensional array with dimensions num_predictors x num_basis x num_basis, where num_predictors is the number of functional predictors, and num_basis is the number of basis functions used to represent each predictor.
#'
#' @param basis_objs A list of basis function objects for each functional predictor
#' @param num_predictors The number of functional predictors
#' @param num_basis The number of basis functions used to represent each predictor
#'
#' @return A three-dimensional array with dimensions num_predictors x num_basis x num_basis
#'
#' @examples
#' # Generate some example data
#' library(fda)
#' data("growth")
#' basis <- create.bspline.basis(rangeval = range(growth$age), nbasis = 10)
#' fdParobj <- fdPar(basisobj = basis, Lfdobj = 2)
#' fd <- smooth.basis(growth$age, growth$hgt, fdParobj)
#' 
#' # Compute the J matrix
#' basis_objs <- list(fd)
#' num_predictors <- 1
#' num_basis <- 10
#' J <- compute_J_matrix_generic(basis_objs, num_predictors, num_basis)
#' 
compute_J_matrix_generic <- function(basis_objs, num_predictors, num_basis) {
  # Initialize the J matrix
  J <-  array(0, dim = c(num_predictors,num_basis,num_basis))

  # Loop over each functional predictor
  for (m in 1:num_predictors) {
    # insert the inner products into the J matrix
    inner <- fda::inprod(basis_objs[[m]], basis_objs[[m]])
    # insert the inner products into the J matrix for the mth predictor
    J[m,,] <- inner
  }
  # Return the J matrix
  return(J)
}




#' Compute Z matrix
#'
#' This function computes the Z matrix used in functional data regression models.
#' The Z matrix is computed as the Kronecker product of the W matrix and the J matrix,
#' where W is an array of predictor functions and J is an array of basis functions.
#' 
#' @param W_array An array of predictor functions.
#' @param J An array of basis functions.
#' @param predictors The number of predictors.
#' @param basis_functions The number of basis functions.
#' 
#' @return A 3-dimensional array representing the Z matrix.
#' 
#' @examples
#' # Define W and J arrays
#' W_array <- array(rnorm(12), dim = c(3, 2, 2))
#' J <- array(rnorm(8), dim = c(2, 2, 2))
#' 
#' # Compute Z matrix
#' Z_matrix <- compute_Z_matrix_generic(W_array, 2, 2, 2)
#' 
#' # Print Z matrix
#' print(Z_matrix)
#'
compute_Z_matrix_generic <- function(W_array, J, predictors, basis_functions) {
  observations <- dim(W_array)[1]
  
  # Initializing Z_matrix
  Z_matrix <- matrix(0, nrow = observations, ncol = predictors * basis_functions)

  # Looping through each observation
  for (i in 1:observations) {
    # Counter to keep track of the column index in Z_matrix
    col_index <- 1 
    for (m in seq_len(predictors)) {
      # Compute the product of W and J for each predictor
      wj_i <- W_array[i,m , ] %*% J[m, , ] 
      # Insert the product into the Z_matrix
      Z_matrix[i, col_index:(col_index + basis_functions - 1)] <- wj_i
      # Update the column index
      col_index <- col_index + basis_functions
    }
  }
  Z_matrix <- array(Z_matrix, dim = c(observations,predictors,basis_functions))
  return(Z_matrix)
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
  RY <- numeric(observations)
  Epsilon <- numeric(observations)
  Ryi <- max(Y) - min(Y)  # Calculate the range of Y
  for (i in 1:observations) {
    RY[i] <- Ryi
    Epsilon[i] <- rnorm(1, mean = 0, sd = sqrt((error_sd * Ryi)^2))
  }
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
compute_noise <- function(Y, error_sd = 0.00) {
  observations <- length(Y)
  noise <- rnorm(observations, mean = 0, sd = error_sd)
  return(noise)
}