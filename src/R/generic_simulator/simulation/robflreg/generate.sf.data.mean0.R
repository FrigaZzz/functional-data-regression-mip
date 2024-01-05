library(robflreg)
library(fda.usc)
# https://github.com/cran/robflreg/blob/master/DESCRIPTION
generate.sf.data.mean0 <-
  function(n, n.pred, n.gp, out.p = 0){
    
    if(!n.pred > 0)
      stop("Error!! Number of predictors must be greater than one!")
    if(!n > 2)
      stop("Error!! Functional variables must have at least two observations!")
    if(!n.gp > 4)
      stop("Error!! The length of grid points must be at least 4 !")
    if(out.p < 0 | out.p > 1)
      stop("Error!! Outlier percentage must be between 0 and 1!")
    
    
    gpX <- seq(0, 1, length.out = n.gp)
    
    cX <- runif(1, min = 1, max = 4)
    
    fX <- fXd <- list()
    for(j in 1:n.pred){
      
      ksi <- list()
      for(i in 1:5){
        ksi[[i]] <- rnorm(n, 1, sd = (cX*i^(-3/2)))
      }
      
      phi <- list()
      for(i in 1:5){
        phi[[i]] <- sin(i * pi * gpX) - cos(i * pi * gpX)
      }
      
      fX[[j]] <- Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
      fXd[[j]] <- Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
    }
    
    coef.space <- list()
    coef.space[[1]] <- sin(pi* gpX)
    coef.space[[2]] <- sin(2*pi* gpX)
    coef.space[[3]] <- sin(3*pi* gpX)
    coef.space[[4]] <- sin(4*pi* gpX)
    coef.space[[5]] <- sin(5*pi* gpX)
    coef.space[[6]] <- cos(pi* gpX)
    coef.space[[7]] <- cos(2*pi* gpX)
    coef.space[[8]] <- cos(3*pi* gpX)
    coef.space[[9]] <- cos(4*pi* gpX)
    coef.space[[10]] <- cos(5*pi* gpX)
    
    coef.ind <- numeric()
    if(n.pred <= 10){
      coef.ind <- sample(1:10, n.pred, replace = FALSE)
    }else{
      coef.ind <- sample(1:10, n.pred, replace = TRUE)
    }
    
    vBeta <- coef.space[coef.ind]
    vBeta0 <- vBeta
    
    for(ij in 1:n.pred){
      fX[[ij]] = fdata(fX[[ij]], argvals = gpX)
      vBeta[[ij]] = fdata(runif(1, min = 1, max = 3) * 
                            vBeta[[ij]], argvals = gpX)
    }
    
    err = rnorm(n, mean=0, sd=1)
    
    fY = Reduce("+", lapply(1:length(fX), function(k){inprod.fdata(fX[[k]], vBeta[[k]])}))
    fYe = fY + err
    
    out.indx <- NULL
    
    if(out.p > 0){
      fX.out <- fXd.out <- list()
      for(j in 1:n.pred){
        
        ksi.out <- list()
        for(i in 1:5){
          ksi.out[[i]] <- rnorm(n, 1, sd = (cX*i^(-1/2)))
        }
        
        phi.out <- list()
        for(i in 1:5){
          phi.out[[i]] <- 2*sin(i * pi * gpX) - cos(i * pi * gpX)
        }
        
        fX.out[[j]] <- Reduce("+", lapply(1:5, function(k){ksi.out[[k]] %*% t(phi.out[[k]])}))
        fXd.out[[j]] <- Reduce("+", lapply(1:5, function(k){ksi.out[[k]] %*% t(phi.out[[k]])}))
      }
      
      coef.ind.out <- numeric()
      if(n.pred <= 10){
        coef.ind.out <- sample(1:10, n.pred, replace = FALSE)
      }else{
        coef.ind.out <- sample(1:10, n.pred, replace = TRUE)
      }
      
      vBeta.out <- coef.space[coef.ind.out]

      for(ij in 1:n.pred){
        fX.out[[ij]] = fdata(fXd.out[[ij]], argvals = gpX)
        vBeta.out[[ij]] = fdata(runif(1, min = 3, max = 5) * 
                                  vBeta.out[[ij]], argvals = gpX)
      }
      
      err = rnorm(n, mean=0, sd=1)
      
      fY.out = Reduce("+", lapply(1:length(fX), function(k){inprod.fdata(fX.out[[k]], vBeta.out[[k]])}))
      fYe.out = fY.out + err
      
      nout <- round(n * out.p)
      out.indx <- sample(1:n, nout)
      
      fYe[out.indx,] <- fYe.out[out.indx,]
      for(io in 1:n.pred)
        fXd[[io]][out.indx,] <- fXd.out[[io]][out.indx,]
    }
      
    
    return(list("Y" = fYe, "X" = fXd, "f.coef" = vBeta0, out.indx = out.indx))
  }

  sim.data <- generate.sf.data.mean0(n = 5, n.pred = 3, n.gp = 5)
