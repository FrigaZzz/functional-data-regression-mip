
measurements <- 50
predictors <- 6
# Define time domains for each ui function as specified in the paper
time_domains <- list(
  seq(0, 1, length.out = measurements),
  seq(0, pi / 3, length.out = measurements), 
  seq(-1, 1, length.out = measurements), 
  seq(0, pi / 3, length.out = measurements), 
  seq(-2, 1, length.out = measurements), 
  seq(-1, 1, length.out = measurements) 
)

# Create the random variables needed for the ui functions as per the paper's description
random_vars <- list(
  a1 = rnorm(1, mean = -4,  sd = sqrt(3 ^ 2)),
  a2 = rnorm(1, mean = 7, sd = sqrt(1.5 ^ 2)),
  b1 = runif(1, min = 3, max = 7),
  b2 = rnorm(1, mean = 0, sd = 1),
  c1 = rnorm(1, mean = -3, sd = sqrt(1.2^2)), 
  c2 = rnorm(1, mean = 2, sd = sqrt(0.5^2)), 
  c3 = rnorm(1, mean = -2, sd = 1), 
  d1 = rnorm(1, mean = -2, sd = 1),
  d2 = rnorm(1, mean = 3, sd = sqrt(1.5^2)),
  e1 = runif(1, min = 2, max = 7),
  e2 = rnorm(1, mean = 2, sd = sqrt(0.4^2)), 
  f1 = rnorm(1, mean = 4, sd = sqrt(2^2)), 
  f2 = rnorm(1, mean = -3, sd = sqrt(0.5^2)), 
  f3 = rnorm(1, mean = 1, sd = 1)
)
u_funcs = list(
  function(t) cos(2 * pi * (t - random_vars$a1)) + random_vars$a2,
  function(t) random_vars$b1 * sin(pi * t) + random_vars$b2,
  function(t) random_vars$c1 * t^3 + random_vars$c2 * t^2 + random_vars$c3 * t,
  function(t) sin(2 * (t - random_vars$d1)) + random_vars$d2 * t,
  function(t) random_vars$e1 * cos(2 * t) + random_vars$e2 * t,
  function(t) random_vars$f1 * exp(-t / 3) + random_vars$f2 * t + random_vars$f3
)

c = 0 # 0.4 or 0.8
beta_funcs = list(
  function(t) {
    sin(t)
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    -c * t^2
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    c * sin(pi * t)
  },
  function(t) {
    0 *t # This function always returns 0 regardless of the input t
  }
)
# Assuming you've defined all other necessary variables and functions
# Now you can run the simulation
set.seed(10)  # for reproducibility
c <- 0  # This should be changed as per the different scenarios you are testing (0, 0.4, or 0.8)
observations <- 100  # Or whatever your number of observations needs to be
predictors <- length(time_domains)  # Based on the number of time domains you have defined



# Placeholder for the final responses
Y <- numeric(observations)


# Initialize the matrices and vectors
FX <- array(NA, dim = c(observations, predictors, measurements))
GU <- matrix(NA, nrow = observations, ncol = predictors)
RY <- numeric(observations)
Epsilon <- matrix(NA, nrow = observations, ncol = 1)
Y <- numeric(observations)


# Fill FX with functional covariate data and calculate RY and Epsilon
for (i in 1:observations) {
  for (m in 1:predictors) {
    time_domain <- time_domains[[m]]
    u_values <- u_funcs[[m]](time_domain)
    RY[i] <- max(u_values) - min(u_values)  # Calculate range for each observation and predictor
    Epsilon[i, ] <- rnorm(1, mean = 0, sd = sqrt(0.025 * RY[i])^2)  # Generate the error term
    FX[i, m, ] <- u_values + Epsilon[i, ]  # Add the error term to the functional covariate

  }
}

# Compute GU by integrating FX * Beta over time_domain
for (i in 1:observations) {
  for (m in 1:predictors) {
    time_domain <- time_domains[[m]]
    beta_func <- beta_funcs[[m]]
    integrand <- function(t) splinefun(time_domain, FX[i, m, ])(t) * beta_func(t)
    GU[i, m] <- integrate(integrand, min(time_domain), max(time_domain))$value
  }
}

# Compute RY and the final response Y
for (i in 1:observations) {
  Ryi <- max(GU[i, ]) - min(GU[i, ])
  RY[i] <- Ryi
  Epsilon[i, ] <- rnorm(1, mean = 0, sd = sqrt((0.05 * Ryi)^2))
  Y[i] <- sum(GU[i, ]) + Epsilon[i, ]
}

