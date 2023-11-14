library(MASS)


# Simulate predictor data from paper
simulate_data_paper <- function(observations, predictors, time_domains, u_funcs) {
  RY <- numeric(observations)
  Epsilon <- matrix(NA, nrow = observations, ncol = 1)
  FX <- array(NA, dim = c(observations, predictors, measurements))

  for (i in 1:observations) {
    for (m in 1:predictors) {
      time_domain <- time_domains[[m]]
      u_values <- u_funcs[[m]](time_domain)
      RY[i] <- max(u_values) - min(u_values)  # Calculate range for each observation and predictor
      Epsilon[i] <- rnorm(1, mean = 0, sd = sqrt(0.025 * RY[i])^2)  # Generate the error term
      FX[i, m, ] <- u_values #+ Epsilon[i]  # Add the error term to the functional covariate
    }
  }
  
  return(list(RY = RY, Epsilon = Epsilon, FX = FX))
}

