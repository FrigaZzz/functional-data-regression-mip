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
  Y <- rep(beta_0, observations)

  # Compute the sum of products of Z_matrix elements and B_matrix coefficients
  for (i in 1:observations) {
    Y[i] <-  sum(Z_matrix[i,,] * B_matrix)
  }

  return(Y)
}

# Compute W matrix, Different time domains
compute_W_matrix_generic <- function(FX,num_basis, time_domains, basis_obj_list) {
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



# Compute the J matrix, Different time domains
compute_J_matrix_generic <- function(basis_objs,num_predictors,num_basis){
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




# Compute Z matrix
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


compute_noise <- function(Y, error_sd = 0.00) {
  observations <- length(Y)
  noise <- rnorm(observations, mean = 0, sd = error_sd)
  return(noise)
}