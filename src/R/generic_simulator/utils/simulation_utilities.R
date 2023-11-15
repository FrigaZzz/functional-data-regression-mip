library(MASS)

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



# Simulate data
simulate_data <- function(mu_func, cov_func_data, time_domain, n) {
  mu <- mu_func(time_domain)
  # extract sig2, rho, decay_type from cov_func_data
  sig2 <- cov_func_data[[1]]
  rho <- cov_func_data[[2]]
  decay_type <- cov_func_data[[3]]
  
  Sigma <- generate_covariance_function(sig2, rho, decay_type)(time_domain, time_domain)
  MASS::mvrnorm(n, mu, Sigma)
}


# Simulate functional features
simulate_functional_features <- function(mu_funcs, cov_funcs_data, n, time_domains) {
  X_list = Map(simulate_data, mu_funcs, cov_funcs_data, time_domains, MoreArgs = list(n = n))

  predictors = length(mu_funcs)
  X <- array(0, dim = c(n, predictors,  length(time_domains[[1]])))

  # Fill the 3D array with data from simulated_features
  for (i in 1:predictors) {
    X[, i, ] <- X_list[[i]]
  }
  return(X)
}


# Simulate predictor data from paper
apply_amplitude_norm <- function(X, observations, predictors) {
  RY <- numeric(observations)
  Epsilon <- matrix(NA, nrow = observations, ncol = 1)
  FX <- array(NA, dim = c(observations, predictors, measurements))

  for (i in 1:observations) {
    for (m in 1:predictors) {
      u_values <- X[i, m, ]
      RY[i] <- max(u_values) - min(u_values)  # Calculate range for each observation and predictor
      Epsilon[i] <- rnorm(1, mean = 0, sd = sqrt(0.025 * RY[i])^2)  # Generate the error term
      FX[i, m, ] <- u_values + Epsilon[i]  # Add the error term to the functional covariate
    }
  }
  
  return(FX)
}


