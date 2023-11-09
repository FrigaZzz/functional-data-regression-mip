library(MASS)

# Simulate data
simulate_data <- function(mu_func, cov_func, n, time) {
  mu <- mu_func(time)
  Sigma <- cov_func(time, time)
  MASS::mvrnorm(n, mu, Sigma)
}


# Simulate response
simulate_response <- function(X, beta_func, cov_func, n, measurements, time) {
  beta_values <- beta_func(time)
  response <- MASS::mvrnorm(n, rep(0, measurements), cov_func(time, time)) + X %*% t(beta_values) / measurements
  return(response)
}


# Simulate functional features
simulate_functional_features <- function(mu_funcs, cov_funcs, n, time) {
  Map(simulate_data, mu_funcs, cov_funcs, MoreArgs = list(n = n, time = time))
}

# Simulate responses for each functional feature
simulate_responses <- function(X_list, beta_funcs, cov_funcs, n, measurements, time) {
  Map(simulate_response, X_list, beta_funcs, MoreArgs = list(cov_func = cov_funcs[[1]], n = n, measurements = measurements, time = time))
}

# Calculate overall mean of Z and binary outcomes
calculate_binary_outcomes <- function(Z_list, mean_Z) {
  sapply(Z_list, function(Z) ifelse(rowMeans(Z) > mean_Z, 1, 0))
}


