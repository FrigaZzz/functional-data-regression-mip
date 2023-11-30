annualprec <- log10(apply(CanadianWeather$dailyAv[,,"Precipitation.mm"],
                          2,sum))
# The simplest 'fRegress' call is singular with more bases
# than observations, so we use a small basis for this example
n_days =  365 #measurements
n_stations = 35
time_points <- seq(0, n_days, length.out = n_days)  # Domain for t
n_basis = 25
smallbasis  <- create.bspline.basis(c(0,n_days), n_basis)


# There are other ways to handle this,
# Simulate temperature data with a U-shape and added variability
temp_data <- matrix(nrow = n_days, ncol = n_stations)

# Assuming the U-shape is symmetric and centered, we can set the center at the midpoint
center <- n_days / 2

for (i in 1:n_stations) {
  # Adding variability to the width of the U-shape
  width <- runif(1, min = 100, max = 150)
  
  # Create a U-shaped curve with added random noise
  temp_data[, i] <- - ((time_points - center) / width)^2 
  temp_data[, i] <- temp_data[, i] + rnorm(n_days, mean = 0, sd = 2)  # Add noise
}

# Convert to functional data
tempfd <- smooth.basis(time_points, temp_data, smallbasis)$fd

# Plot the functional data
plot(tempfd, xlab = "Day", ylab = "Temperature")



xfdlist <- list(const=rep(1, n_stations), tempfd=tempfd)

# The intercept must be constant for a scalar response
betabasis1 <- create.constant.basis(c(0, n_days))
betafd1    <- fd(0, betabasis1)
betafdPar1 <- fdPar(betafd1)

betafd2     <- with(tempfd, fd(basisobj=smallbasis, fdnames=fdnames))
# convert to an fdPar object
betafdPar2  <- fdPar(betafd2)

betalist <- list(const=betafdPar1, tempfd=betafdPar2)

precip_out<- fRegress(annualprec, xfdlist, betalist)
