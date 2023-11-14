library(fda)


# Compute W matrix
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

# Compute W matrix, Different time domains
compute_W_matrix_paper <- function(FX,num_basis, time_domains, basis_obj_list) {
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


# Compute the J matrix (inner products for basis functions)
compute_J_matrix <- function(basis_obj) {
  fda::inprod(basis_obj, basis_obj)
}


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





# # Compute Y matrix
compute_Y_values <- function(observations, predictors, basis_functions, B_matrix, beta_0, Z_matrix) {
  Y <- numeric(observations)
  for (i in 1:observations) {
    Y[i] <- beta_0  # Adding the intercept first
    for (m in 1:predictors) {
      for (j in 1:basis_functions) {
        Y[i] <- Y[i] + B_matrix[m, j] * Z_matrix[i,m,j]
      }
    }
  }
  return(Y)
}


compute_noise <- function(Y, error_sd) {
  observations <- length(Y)
  noise <- rnorm(observations, mean = 0, sd = error_sd)
  return(noise)
}

