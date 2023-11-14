library(MASS)

# Simulate data
simulate_data <- function(mu_func, cov_func, n, time) {
  mu <- mu_func(time)
  Sigma <- cov_func(time, time)
  MASS::mvrnorm(n, mu, Sigma)

}


# Simulate functional features
simulate_functional_features <- function(mu_funcs, cov_funcs, n, time) {
  X_list = Map(simulate_data, mu_funcs, cov_funcs, MoreArgs = list(n = n, time = time))

  predictors = length(mu_funcs)
  X <- array(0, dim = c(n, predictors,  length(time)))

  # Fill the 3D array with data from simulated_features
  for (i in 1:predictors) {
    X[, i, ] <- X_list[[i]]
  }
  return(X)
}


# Simulate predictor data from paper
simulate_data_paper <- function(observations, predictors, time_domains, u_funcs) {
  RY <- numeric(observations)
  Epsilon <- matrix(NA, nrow = observations, ncol = 1)
  FX <- array(NA, dim = c(observations, predictors, measurements))

  for (i in 1:observations) {
    for (m in 1:predictors) {
      time_domain <- time_domains[[m]]
      u_values <- u_funcs[[m]](time_domain)
      FX[i, m, ] <- u_values 
    }
  }
  return(FX = FX)
}


# Simulate predictor data from paper
apply_amplitude_norm <- function(FX, observations, predictors) {
  RY <- numeric(observations)
  Epsilon <- matrix(NA, nrow = observations, ncol = 1)
  FX <- array(NA, dim = c(observations, predictors, measurements))

  for (i in 1:observations) {
    for (m in 1:predictors) {
      u_values <- FX[i, m, ]
      RY[i] <- max(u_values) - min(u_values)  # Calculate range for each observation and predictor
      Epsilon[i] <- rnorm(1, mean = 0, sd = sqrt(0.025 * RY[i])^2)  # Generate the error term
      FX[i, m, ] <- u_values + Epsilon[i]  # Add the error term to the functional covariate
    }
  }
  
  return(FX = FX)
}


