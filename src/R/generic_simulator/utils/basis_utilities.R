# Load necessary libraries
library(fda)

create_1_basis <- function(basis_functions, bm_range, norder) {
  # Loop over each functional predictor
    # Use the specific time domain for each functional predictor to create  the basis object
    basis_obj<- fda::create.bspline.basis(rangeval = c(min (bm_range), max(bm_range)), nbasis = basis_functions, norder = norder)
  return(basis_obj)
}
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

# Function to expand functional data into B-spline basis
expand_functional_data <- function(data_matrix, time, basis_obj) {
  fd_list <- apply(data_matrix, 1, function(row) fda::smooth.basis(time, row, basis_obj)$fd)
  do.call(rbind, lapply(fd_list, function(fd) t(fd$coefs)))
}

# Smooth beta coefficients into the B-spline basis
smooth_beta_function <- function(beta_func,time,basis_obj) {
  beta_values <- beta_func(time)
  fdPar_obj <- fdPar(basis_obj)
  smoothed_beta <- smooth.basis(time, beta_values, fdPar_obj)
  beta_coefs <- coef(smoothed_beta$fd)
  beta_coefs <- beta_coefs[, 1]
  return(beta_coefs)
}

smooth_betas <- function(beta_funcs,num_basis, time, basis_obj) {
  num_predictors <- length(beta_funcs)
  
  # Initialize a matrix to store the smoothed beta coefficients
  beta_matrix <-  array(0, dim = c(num_predictors,num_basis))
  beta_point_values <-  array(0, dim = c(num_predictors,length(time)))
  basis_values <-   eval.basis(time, basis_obj)
  for (i in 1:num_predictors) {
    beta_func <- beta_funcs[[i]]    
    beta_values <- beta_func(time)
    beta_point_values[i,] <- beta_values
    fdPar_obj <- fda::fdPar(basis_obj)
    smoothed_beta <- fda::smooth.basis(time, beta_values, fdPar_obj)
    beta_matrix[i,] <- smoothed_beta$fd$coefs
  }
  return(list(beta_matrix = beta_matrix, beta_point_values =   beta_point_values, basis_values = basis_values))
}

# Smooth beta coefficients into the B-spline basis, Different time domains
smooth_betas_generic <- function(beta_funcs,num_basis, time_domains, basis_objs) {
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
