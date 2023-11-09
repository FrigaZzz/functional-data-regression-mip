# Load necessary libraries
library(fda)

# Function to expand functional data into B-spline basis
expand_functional_data <- function(data_matrix, time, basis_obj) {
  fd_list <- apply(data_matrix, 1, function(row) fda::smooth.basis(time, row, basis_obj)$fd)
  do.call(rbind, lapply(fd_list, function(fd) fd$coefs))
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
