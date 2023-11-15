# Covariance function generator

#' Generate a covariance function for a given set of parameters
#'
#' @param sig2 The variance parameter
#' @param rho The range parameter
#' @param decay_type The type of decay function to use (either "exponential" or "matern")
#'
#' @return A covariance function that takes two arguments (t and s) and returns a covariance matrix
#'
generate_covariance_function <- function(sig2, rho, decay_type = "exponential") {
  function(t, s) {
    d <- abs(outer(t, s, "-"))
    if (decay_type == "exponential") {
      return(sig2 * exp(-d / rho))
    } else if (decay_type == "matern") {
      return(sig2 * (1 + sqrt(5) * d / rho) * exp(-sqrt(5) * d / rho))
    }
  }
}


test = FALSE
if (test) {
  # Generate a covariance function with exponential decay
  cov_func_exp <- generate_covariance_function(sig2 = 1, rho = 0.5, decay_type = "exponential")

  # Generate a covariance function with Matérn decay
  cov_func_matern <- generate_covariance_function(sig2 = 1, rho = 0.5, decay_type = "matern")

  # Define some time points
  t <- seq(0, 1, length.out = 100)
  s <- seq(0, 1, length.out = 100)

  # Compute the covariance matrices
  cov_matrix_exp <- cov_func_exp(t, s)
  cov_matrix_matern <- cov_func_matern(t, s)

  # Print the covariance matrices
  print(cov_matrix_exp)
  print(cov_matrix_matern)
}