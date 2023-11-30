library(fda) 
set.seed(123)

n_observations <- 35
n_measurements <- 365  # Assuming daily measurements over a year
n_basis <- 6

# Basis for functional data
basis <- create.fourier.basis(c(0, 365), n_basis)

# Coefficients for predictors
a1 <- rnorm(n_observations, mean=-4, sd=3)
a2 <- rnorm(n_observations, mean=1.5, sd=0.5)
b1 <- runif(n_observations, min=3, max=7)
b2 <- rnorm(n_observations, mean=0, sd=1)

# Time points
t <- seq(0, pi/3, length.out = n_measurements)

# Predictor functions
u1 <- matrix(0, n_measurements, n_observations)
u2 <- matrix(0, n_measurements, n_observations)

for (i in 1:n_observations) {
  u1[, i] <- a1[i] * cos(2 * pi * (t / 365)) + a2[i]
  u2[, i] <- b1[i] * sin(2 * pi * (t / 365)) + b2[i]
}

# Smooth the predictors
fd1 <- smooth.basis(t, u1, basis)$fd
fd2 <- smooth.basis(t, u2, basis)$fd

plot(fd1)

# Create fdPar objects for regression
beta_fdPar1 <- fdPar(fd1, Lfdobj=int2Lfd(2), lambda=0)
beta_fdPar2 <- fdPar(fd2, Lfdobj=int2Lfd(2), lambda=0)

# List of fdPar objects for regression
betalist <- list(beta_fdPar1, beta_fdPar2)

# Annual precipitation (response variable)
# Here we sum the values of u1 as a proxy for precipitation
annualprec <- apply(u1, 2, sum)

# Perform functional regression
precip <- fRegress(annualprec ~ fd1 + fd2)
