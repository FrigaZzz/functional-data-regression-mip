library(fda)

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
  num_observations <- dim(FX)[1]
  num_predictors <- dim(FX)[2]
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
      wj_i <- t(W_array[i,m , ]) %*% J[m, , ] 
      # Insert the product into the Z_matrix
      Z_matrix[i, col_index:(col_index + basis_functions - 1)] <- wj_i
      # Update the column index
      col_index <- col_index + basis_functions
    }
  }
  Z_matrix_out <- array(0, dim = c(observations,predictors,basis_functions))
  # fill Z matrix 
  for (i in 1:observations) {
    for (m in 1:predictors) {
      Z_matrix_out[i,m,] <- Z_matrix[i,((m-1)*basis_functions+1):(m*basis_functions)]
    }
  }
  return(Z_matrix_out)
}


