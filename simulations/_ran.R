require(splines)
require(grplasso)

grplFlinear <- function(Y, X, Tps, lambda, phi, dfs = 10,
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
  
  print(dim(Y))
  ## group lasso
  index = c(NA,rep(1:nfunc,dfs))
  grpl = grplasso(x = Z, y = Y, index = index, model = LinReg(), lambda = lambda,
                  standardize = F, ...)
  
  
  ## output: intercept and fitted coefficient functions
  intercept = grpl$coef[1,]
  Coef <- list()
  index[1] = 0
  for (jj in 1:nfunc)
  {
    Coef[[jj]] <- B[[jj]]%*%iR[[jj]]%*%grpl$coef[index == jj,]
  }
  
  out = list("intercept" = intercept, "Coef" = Coef)
  return(out)
}

simfx1 <- function(n, p, tps, func_shapes, varx = rep(0, p), coef_functions, noise_sd = 1, exclude_predictor = NULL) {
  fx <- vector("list", p)
  fobs <- vector("list", p)
  
  simulate_predictor <- function(n, tps, func_shape, varx, bij, mij) {
    fx_single <- matrix(0, n, length(tps))
    fobs_single <- matrix(0, n, length(tps))
    
    for (i in 1:n) {
      fx_single[i, ] <- sapply(tps, func_shape, bij = bij, mij = mij)
      fobs_single[i, ] <- fx_single[i, ] + rnorm(length(tps), 0, sqrt(varx))
    }
    
    list("fx" = fx_single, "fobs" = fobs_single)
  }
  
  if (!is.list(func_shapes) || length(func_shapes) != p || 
      !is.list(coef_functions) || length(coef_functions) != p) {
    stop("func_shapes and coef_functions must be lists with length equal to p")
  }
  
  Y <- rep(0, n)
  for (j in 1:p) {
    # Generate coefficients from a normal distribution
    bij <- rnorm(1, mean = 0, sd = 1)  # Adjust mean and sd as needed
    mij <- rnorm(1, mean = 0, sd = 1)  # Adjust mean and sd as needed
    
    result <- simulate_predictor(n, tps[[j]], func_shapes[[j]], varx[j], bij, mij)
    fx[[j]] <- scale(result$fx) / 10
    fobs[[j]] <- result$fobs
    if (is.null(exclude_predictor) || !j %in% exclude_predictor) {
      Y <- Y + rowSums(fobs[[j]] * sapply(tps[[j]], coef_functions[[j]]))
    }
  }
  
  Y <- Y + rnorm(n, 0, noise_sd)
  
  list("funx" = fx, "funcs" = fobs, "Y" = Y)
}

# Example usage remains the same

# Example Usage
set.seed(123)
tps <- list(1:300, 1:300, 1:300)
tp <- 1:300
n <- 400
func_shapes <- list(
  function(tp, bij, mij) (sum(bij * sin(tp * (5 - bij) * (2 * pi / max(tps[[1]]))) - mij) + 15) / 100,
  function(tp, bij, mij) (sum(bij * cos(tp * (5 - bij) * (2 * pi / max(tps[[1]]))) - mij) + 15) / 100,
  function(tp, bij, mij) (sum(bij * tp*(tp * (5 - bij) * (2 * pi / max(tps[[1]]))) - mij) + 15) / 100
)
coef_functions <- list(
  function(tp) sin(2 * pi * tp / max(tps[[1]])),
  function(tp) 3*(tp+1),
  function(tp) rep(0, length(tp))
)
output <- simfx1(n = n, p = 3, tps = tps, func_shapes = func_shapes, 
                 varx = rep(0.1, 3), coef_functions = coef_functions, 
                 noise_sd = 2, exclude_predictor = c(3))

fx  = output$funx
fobs = output$funcs
Y = output$Y

lambda <- 10^seq(3,0,by=-1)
phi <- 10^c(10,8,6,4,2)
grpl1 <- grplFlinear(Y = Y, X = fobs, Tps = tps, lambda = lambda, phi = phi[1],
                     dfs = 6)

# Assuming you have already run the grplFlinear and simfx1 functions

# True Coefficients
true_coef1 <- sapply(tps[[1]], coef_functions[[1]])
true_coef2 <- sapply(tps[[2]], coef_functions[[2]])
true_coef3 <- sapply(tps[[2]], coef_functions[[3]])

# Estimated Coefficients
est_coef1 <- grpl1$Coef[[1]][,1]
est_coef2 <- grpl1$Coef[[2]][,1]
est_coef3 <- grpl1$Coef[[3]][,1]




plot(tp, rep(0,length(tp)), type="l", ylim=c(-0.2,0.2), xlab="t", ylab=expression(beta[2](t)), bty="n")
for(wl in 1:length(lambda))
  lines(tp,grpl1$Coef[[2]][,wl],lty=wl)
title(expression(varphi==10^10))
