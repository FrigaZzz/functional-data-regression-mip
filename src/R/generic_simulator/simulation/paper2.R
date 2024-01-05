library(MASS)
library(here)
library(fda)

simulate_paper2_data <- function(observations, time_domains, predictors){
  tps <- time_domains
  n <- observations
  XX <- simfx1(n = n, p = predictors, tps = tps, varx = rep(0,predictors))

  # true coefficient functions
  measurements <- length(tps[[1]])
  tp <- 1:measurements
  bt <- 5*dgamma(tp/10, 3, 1/3)
  # y values
  sigma2 <- 4
  mu <- XX$funx[[1]]%*%bt
  Y <- mu + rnorm(n, 0, sqrt(sigma2))

  # reshape X and U
   X <- array(dim = c(observations, predictors, measurements))
   U <- array(dim = c(observations, predictors, measurements))

   for (j in 1:predictors) {
     X[, j, ] <- (XX$funx[[j]])
     U[, j, ] <- (XX$funcs[[j]])
   }

  # X = XX$funx
  # U = XX$funcs

  return(list("X" = X, "U" = U, "Y" = Y))
}

simfx1 <- function(n, p, tps, varx=rep(0,p), bx=5, mx=2*pi) {
  fx <- list()
  fobs <- list()
  # all_bij <- lapply(1:p, function(x) matrix(runif(5 * n, 0, bx), n, 5))
  # all_mij <- lapply(1:p, function(x) matrix(runif(5 * n, 0, mx), n, 5))
  for (j in 1:p) {
    tmax <- max(tps[[j]])
    fx[[j]] <- matrix(0, n, length(tps[[j]]))
    fobs[[j]] <- matrix(0, n, length(tps[[j]]))
    for (i in 1:n) {
      bij <- runif(5, 0, bx)
      mij <- runif(5, 0, mx)

      tfx <- function(tp) {
        (sum(bij*sin(tp*(5-bij)*(2*pi/tmax)) - mij) + 15)/100
      }
      fx[[j]][i,] <- sapply(tps[[j]], tfx)
    }
  }
  fx <- lapply(fx, scale)
  for (j in 1:p) {
    fx[[j]] <- fx[[j]]/10
    for (i in 1:n) {
      fobs[[j]][i,] <- fx[[j]][i,] + rnorm(length(tps[[j]]), 0, sqrt(varx[j]))
    }
  }
  return(list("funx"=fx, "funcs"=fobs))
}


# set.seed(1)

# # run test 
# time_domains <- list(
#   list(0, 300),
#   list(0, 300)

# )

# time_domains_eval <- lapply(time_domains, function(domain) {
#             seq(from = domain[[1]], to = domain[[2]], length.out = 300)
#         })


# y <- simulate_paper2_data(observations = 300, time_domains = time_domains_eval, predictors = 2)

# # plot(y$X[[1]][1,])