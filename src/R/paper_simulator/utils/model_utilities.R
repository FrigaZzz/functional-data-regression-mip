library(fda)

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



# Compute the J matrix, Different time domains
compute_J_matrix_paper <- function(basis_objs,num_predictors,num_basis){
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
compute_Z_matrix_paper <- function(W_array, J, predictors, basis_functions) {
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



# # Compute Y matrix
compute_Y_values <- function(observations, predictors, basis_functions, B_matrix, beta_0, Z_matrix) {
  Y <- numeric(observations)
  for (i in 1:observations) {
    Y[i] <- beta_0  # Adding the intercept first
    for (m in 1:predictors) {
      sum_B_Z <- 0  # Summation of B and Z products for each predictor
      for (j in 1:basis_functions) {
        sum_B_Z <- sum_B_Z + B_matrix[m, j] * Z_matrix[i,m,j]
      }
      Y[i] <- Y[i] + sum_B_Z  # Adding the summation to Y
    }
  }
  return(Y)
}



compute_noise_paper<- function(Y) {
  observations <- length(Y)
  RY <- numeric(observations)
  Epsilon <- numeric(observations)
  Ryi <- max(Y) - min(Y)  # Calculate the range of Y
  for (i in 1:observations) {
    RY[i] <- Ryi
    Epsilon[i] <- rnorm(1, mean = 0, sd = sqrt((0.05 * Ryi)^2))
  }
  return(Epsilon)
}

