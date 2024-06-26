# Load necessary libraries
library(refund)
library(MASS)
library(fda)
library(here)

# Source utility files
source(here("src","R", "generic_simulator", "config.R")) # sets the utility path
source(here("src", "R",  "generic_simulator",   "simulation" , "cov.R"))
source(here("src", "R",  "generic_simulator",   "simulation" , "paper.R"))
source(here("src", "R",  "generic_simulator",   "simulation" , "paper2.R"))

require(splines)
require(grplasso)


compute_Z_penalized <- function(Y, X, Tps, lambda, phi = 0.01, dfs = 10,
adapt1 = NULL, adapt2 = NULL, ...){

	#### Observed functions x_{ij}(t) at possibly different support grids t
	## Y = nx1 vector of responses
	## X = list of matrices, each element corresponds to observations of one functional predictor.
	##		Thus, [[jj]] would provide the jj-th observed function with subjects as rows and the support grid as columns
	##		Note that each function (i.e. component of X) can have different support
	## Tps = list of vectors, each element is a suppport grid of the corresponding observed function.
	##		Thus, Tps[[jj]] would give the support grid of the jj-th observed function. This grid is
	##		assumed to be the same for all subjects ii, and equidistant
	## lambda = vector of penalty parameters
	## phi = penalty parameter for smoothing
	## dfs = dfs used for basis expansions of coefficient functions (can be a vector)


	nsub = length(Y) ## number of subjects
	nfunc = length(Tps) ## number of functions per subject



	#### We use bsplines as basis functions for the corrsponding beta functions
  if (length(dfs) == 1)
    dfs = rep(dfs, nfunc) ## vector of intended df of each spline basis
  if (length(dfs) != nfunc)
    stop("length of dfs does not match number of predictors")

	B <- Psi <- Omega <- K <- iR <- eK <- list()
	delt <- rep(NA, nfunc)


	for (jj in 1:nfunc){
		
    spj = diff(range(Tps[[jj]]))#/(dfs[jj]-2)
    bknj = c(min(Tps[[jj]]) - spj, max(Tps[[jj]]) + spj) ## boundary knots
		B[[jj]] = bs(Tps[[jj]], df=dfs[jj], Boundary.knots=bknj) ## basis spline set up
		delt[jj] = Tps[[jj]][2] - Tps[[jj]][1] ## differences in Tps

		Psi[[jj]] = delt[jj] * t(B[[jj]]) %*% B[[jj]] ## approximate norm of bsplines assuming dense design
    if (length(adapt1) == nfunc)
      Psi[[jj]] = adapt1[jj]*Psi[[jj]]

    dBj <- matrix(NA,nrow(B[[jj]]),ncol(B[[jj]]))
    for (k in 1:ncol(B[[jj]])) ## computation of second derivatives
      {
        iS <- interpSpline(Tps[[jj]],B[[jj]][,k])
        dBj[,k] <- predict(iS, Tps[[jj]], deriv = 2)$y
      }
    Omega[[jj]] = delt[jj] * t(dBj) %*% dBj ## approximate norm of 2nd deriv of bsplines assuming dense design
    if (length(adapt2) == nfunc)
      Omega[[jj]] = adapt2[jj]*Omega[[jj]]

		K[[jj]] = Psi[[jj]] + phi * Omega[[jj]] ## K matrix

    # Check eigenvalues
    eig_vals <- eigen(K[[jj]])$values
    if (any(eig_vals <= 0)) {
        cat("Matrix K[", jj, "] has non-positive eigenvalues.\n")
        # Diagonal loading for numerical stability
        K[[jj]] <- K[[jj]] + diag(abs(min(eig_vals)) + 1e-6, ncol(K[[jj]]))
    }


    eK[[jj]] <- eigen(K[[jj]])
    #iR[[jj]] <- t((1/sqrt(eK[[jj]]$values))*t(eK[[jj]]$vectors))
		iR[[jj]] = backsolve(chol(K[[jj]]), x = diag(ncol(K[[jj]])))  ## inverse of cholesky of K
	}


	## covariates for the linear model
	Z = 1
	for (jj in 1:nfunc)
    {
			tmp = delt[jj]*(X[[jj]]%*%B[[jj]])
			Z = cbind(Z, tmp%*%iR[[jj]])
		}
return (list(Z = Z, K = K, iR = iR, eK = eK, B = B))
}


# Unified function
generate_data <- function(
    predictors, observations, measurements, basis_functions, intercept = 0, norder, noise_snr, 
    mu_funcs, beta_funcs, time_domains, cov_funcs = NULL, seed = 2000, phi,
    simulation_type = "not paper") {

  # set.seed(seed)
  # debug print all input parameters
  print(paste("predictors:", predictors))
  print(paste("observations:", observations))
  print(paste("measurements:", measurements))
  print(paste("basis_functions:", basis_functions))
  print(paste("intercept:", intercept))
  

  # Call the appropriate data generation function based on simulation type
  if (simulation_type == "paper") {
    data <- simulate_paper_data(mu_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr)
   } else if (simulation_type == "paper2") {
    data <- simulate_paper2_data( observations, time_domains,  predictors) 
  }else {
    data <- simulate_cov_data(mu_funcs, cov_funcs, beta_funcs, observations, time_domains, intercept, predictors, noise_snr) 
  }

  # Extract U, X, Y from the returned data
  U <- data$U
  X <- data$X
  Y <- data$Y
  beta_point_values <- create_beta_curves(beta_funcs, time_domains)

  fobs <- list()
  for (j in 1:predictors) {
    fobs[[j]] <- matrix(0, observations, measurements)
    for (i in 1:observations) {
      fobs[[j]][i,] <- X[i,j,]
    }
  }
  
  tps <-time_domains
  
  # Remaining processing steps
  basis_objs <- create_basis(basis_functions, time_domains, norder, predictors)
  result <- smooth_betas_generic(beta_point_values, basis_functions, time_domains, basis_objs)
  Beta_matrix <- result$beta_matrix
  basis_values <- result$basis_values
  J <- compute_J_matrix_generic(basis_objs, predictors, basis_functions)
  W <- compute_W_matrix_generic(X, basis_functions, time_domains, basis_objs)
  print(phi )
  z_computation_result = compute_Z_penalized(Y = Y, X = fobs, Tps = tps, phi = phi, dfs = basis_functions )

  Z_matrix = z_computation_result$Z
  K_matrix = z_computation_result$K
  iR_matrix = z_computation_result$iR
  eK_matrix = z_computation_result$eK
  Bs = z_computation_result$B

  list(W = W, Z = Z_matrix, Y = Y, J = J, B = Beta_matrix, U = U, X = X, basis_objs = basis_objs, basis_values = basis_values, beta_point_values = beta_point_values, K = K_matrix, iR = iR_matrix, eK = eK_matrix, Bs = Bs)

}

# # test = TRUE
# source(here("simulations","settings", "3_predictors", "default.R")) # sets the utility path
# time_domains = list(
#   seq(0, 1, length.out = measurements),
#   seq(0, pi / 3, length.out = measurements),
#   seq(-1, 1, length.out = measurements)
# )
# phi = 0.01
# data <- generate_data(predictors, observations,measurements, basis_functions, intercept, norder,noise_snr, mu_funcs, beta_funcs,  time_domains, cov_funcs,seed,phi, simulation_type)

# # Assuming `data` is your data output variable
# W = data$W
# Z = data$Z
# Y = data$Y
# J = data$J
# B = data$B
# U = data$U
# X = data$X
# basis_objs = data$basis_objs
# basis_values = data$basis_values
# beta_dt = data$beta_dt