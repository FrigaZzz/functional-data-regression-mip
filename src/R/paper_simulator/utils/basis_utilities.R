# Load necessary libraries
library(fda)

# Smooth beta coefficients into the B-spline basis, Different time domains
smooth_betas_paper <- function(beta_funcs,num_basis, time_domains, basis_objs) {
  num_predictors <- length(beta_funcs)
  
  # Initialize a matrix to store the smoothed beta coefficients
  beta_matrix <-  array(0, dim = c(num_predictors,num_basis))
  beta_point_values <-  array(0, dim = c(num_predictors,length(time_domains[[1]])))
  basis_values <-  array(0, dim = c(num_predictors,num_basis,length(time_domains[[1]])))
  for (i in 1:num_predictors) {
    beta_func <- beta_funcs[[i]]
    basis_obj <- basis_objs[[i]]
    time <- time_domains[[i]]

    beta_values <- beta_func(time)
    beta_point_values[i,] <- beta_values
    basis_values[i,,] <- eval.basis(time, basis_obj)
    fdPar_obj <- fda::fdPar(basis_obj)
    smoothed_beta <- fda::smooth.basis(time, beta_values, fdPar_obj)
    beta_matrix[i,] <- smoothed_beta$fd$coefs
  }
  return(list(beta_matrix = beta_matrix, beta_point_values = beta_point_values, basis_values = basis_values))

}