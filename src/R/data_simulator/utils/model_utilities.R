library(fda)

# Compute the J matrix (inner products for basis functions)
compute_J_matrix <- function(basis_obj) {
  fda::inprod(basis_obj, basis_obj)
}

# Compute Z matrix
compute_Z_matrix <- function(W_array, J, predictors, basis_functions) {
  observations <- dim(W_array)[1]
  
  # Initializing Z_matrix with a column of ones for the intercept
  Z_matrix <- matrix(0, nrow = observations, ncol = predictors * basis_functions)

  # Looping through each observation
  for (i in 1:observations) {
    # Counter to keep track of the column index in Z_matrix
    col_index <- 1 
    for (m in seq_len(predictors)) {
      for (b in seq_len(basis_functions)) {
        # Calculating and assigning the transformed predictors to the Z_matrix
        Z_matrix[i, col_index] <- W_array[i, m, b] * sum(J[b, ] * W_array[i, m, ])
        col_index <- col_index + 1
      }
    }
  }
  return(Z_matrix)
}

# Compute Y matrix
compute_Y_values <- function(observations, predictors, basis_functions, B_matrix, beta_0, Z_matrix, all_errors) {
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
    Y[i] <- Y[i] + all_errors[i]  # Adding the error term at the end
  }
  return(Y)
}

