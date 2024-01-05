library(robflreg)


simulate_data_robust <- function(n_pred, n_curve, n_gp) {
  # Generate a dataset with n_pred functional predictors and n_curve observations
  sim.data <- generate.sf.data(n = n_curve, n.pred = n_pred, n.gp = n_gp)

  # Extract response variable and predictors
  Y <- sim.data$Y
  X <- sim.data$X
  X <- U <- reshape_X(X)

  # Extract true parameter functions and store them in a 2D array
  beta_values <- array(0, dim = c(n_pred, n_gp))
  for (i in 1:n_pred) {
    beta_values[i, ] <- sim.data$f.coef[[i]]
  }

  # Return a list containing the generated data
  return(list(X = X, U = X, Y = Y, beta_values = beta_values))
}
reshape_X <- function(X) {
  # Get the number of predictors, observations, and measurements
  n_pred <- length(X)
  n_curve <- dim(X[[1]])[1]
  n_gp <- dim(X[[1]])[2]

  # Initialize a 3D array to store the reshaped X
  X_reshaped <- array(0, dim = c(n_curve, n_pred, n_gp))

  # For each predictor
  for (j in 1:n_pred) {
    # Assign the matrix corresponding to the predictor to the appropriate slice of the 3D array
    X_reshaped[, j, ] <- X[[j]]
  }

  # Return the reshaped X
  return(X_reshaped)
}



# set.seed(1)
# # Usage:
# data <- generate_data(n_pred = 5, n_curve = 400, n_gp = 101)
# #plot(data$X[[1]][1,])

library(robflreg)
library(here)
# Source utility files using the relative paths
source(here("src", "R",  "generic_simulator",   "simulation" , "robflreg","generate.sf.data.mean0.R"))

simulate_data_robust <- function(n_pred, n_curve, n_gp) {
  # Generate a dataset with n_pred functional predictors and n_curve observations
  sim.data <- generate.sf.data.mean0(n = n_curve, n.pred = n_pred, n.gp = n_gp)

  # Extract response variable and predictors
  Y <- sim.data$Y
  X <- sim.data$X
  X <- U <- reshape_X(X)

  # Extract true parameter functions and store them in a 2D array
  beta_values <- array(0, dim = c(n_pred, n_gp))
  for (i in 1:n_pred) {
    beta_values[i, ] <- sim.data$f.coef[[i]]
  }

  # Return a list containing the generated data
  return(list(X = X, U = X, Y = Y, beta_values = beta_values))
}
reshape_X <- function(X) {
  # Get the number of predictors, observations, and measurements
  n_pred <- length(X)
  n_curve <- dim(X[[1]])[1]
  n_gp <- dim(X[[1]])[2]

  # Initialize a 3D array to store the reshaped X
  X_reshaped <- array(0, dim = c(n_curve, n_pred, n_gp))

  # For each predictor
  for (j in 1:n_pred) {
    # Assign the matrix corresponding to the predictor to the appropriate slice of the 3D array
    X_reshaped[, j, ] <- X[[j]]
  }

  # Return the reshaped X
  return(X_reshaped)
}



# set.seed(1)
# # Usage:
# data <- generate_data(n_pred = 5, n_curve = 400, n_gp = 101)
# #plot(data$X[[1]][1,])

# # plot(data$beta_values[1,], type = "l", ylim = c(-1, 1))

# recompute_Y <- function(X, beta_values) {
#   # Get the number of predictors, observations, and measurements
#   n_pred <- dim(X)[2]
#   n_curve <- dim(X)[1]
#   n_gp <- dim(X)[3]

#   # Initialize Y_recomputed
#   Y_recomputed <- array(0, dim = c(n_curve))

#   # For each observation
#   for (i in 1:n_curve) {
#     # For each predictor
#     for (j in 1:n_pred) {
#       # Compute the integral of the product of the predictor function and the beta function
#       print(X[i, j, ])
#       Y_recomputed[i] <- Y_recomputed[i] + sum(X[i, j, ] * beta_values[j, ])
#     }
#   }

#   # Return the recomputed Y
#   return(Y_recomputed)
# }

# Y_r= recompute_Y(data$X, data$beta_values)

