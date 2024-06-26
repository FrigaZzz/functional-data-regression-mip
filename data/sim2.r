# ------------------------------------------------------------------------------
# -------- A Simulation Study for grplMFLM -------------------------------------
# -------- dense design, normal y          -------------------------------------
# ------------------------------------------------------------------------------

# data function
library('mgcv')
library('refund')
source("C:\\Users\\User\\Documents\\repos\\TESI\\tesi_data\\progetto_s_on_f.md\\b_spline\\functional-data-regression-mip\\data\\fdatagen.r")
# test
tps <- list()
tps[[1]] <- tps[[2]] <- tps[[3]] <- tps[[4]] <- 1:300
n <- 3
XX <- simfx1(n = n, p = 4, tps = tps, varx = rep(0,4))
# par(mfrow = c(2,2))
# for (j in 1:4)
#   {
#     plot(XX$funcs[[j]][1,], type = "l", ylim = range(XX$funcs[[j]]))
#     for (i in 1:n)
#       lines(XX$funcs[[j]][i,])
#   }
# par(mfrow = c(1,1))

# true coefficient functions
p <- 30
#p <- 10
tp <- 1:300
bt <- dgamma(tp/10, 3, 1/3)
#plot(tp, bt, type="l")
#plot(tp, 5.5*bt, type="l")
btrue <- list()
Tps <- list()
for (j in 1:p)
  {
    btrue[[j]] <- ifelse(j<=6,10-1.5*j,0)*bt
    Tps[[j]] <- tp
  }
btrue[[7]] <- 0.6*dexp(0.02*tp)
btrue[[8]] <- 0.4*dexp(0.02*tp)
btrue[[9]] <- 0.6*dexp(0.01*tp)
btrue[[10]] <- 0.4*dexp(0.01*tp)
  
# expected y
n <- 1000
mu <- rep(10,n)
XX <- simfx1(n = n, p = p, tps = Tps)
for (jj in 1:length(XX$funx))
  {
    mu <- mu + XX$funx[[jj]]%*%btrue[[jj]]
  }
# hist(mu)
var(mu)
kb <- ceiling(n/p) - 1
smoothFit <- pfr_old(Y = mu + rnorm(n,0,2), funcs = XX$funcs, kb = min(30,kb))
# plot(smoothFit$BetaHat[[1]])
# lines(btrue[[1]])
  
  
# grpl function
source("C:\\Users\\User\\Documents\\repos\\TESI\\tesi_data\\progetto_s_on_f.md\\b_spline\\functional-data-regression-mip\\data\\grplMFLM.r")
# simulation
set.seed(111)
n <- 300
#R <- 20
R <- 50

# estimated coefs
intercepts <- numeric(R)
interceptsA1 <- numeric(R)
interceptsA2 <- numeric(R)
interceptsS <- numeric(R)
interceptsW <- numeric(R)

EstCoefs <- list()
EstCoefsAdapt1 <- list()
EstCoefsAdapt2 <- list()
EstCoefsSmooth <- list()
EstCoefsWiggly <- list()
for (j in 1:p)
  {
    EstCoefs[[j]] <- matrix(NA,length(Tps[[j]]),R)
    EstCoefsAdapt1[[j]] <- matrix(NA,length(Tps[[j]]),R)
    EstCoefsAdapt2[[j]] <- matrix(NA,length(Tps[[j]]),R)
    EstCoefsSmooth[[j]] <- matrix(NA,length(Tps[[j]]),R)
    EstCoefsWiggly[[j]] <- matrix(NA,length(Tps[[j]]),R)
}



# error variance etc
sigma2 <- 10
lambda <- 10^seq(3,1,by=-0.1)
phi <- 10^seq(8,4,by=-1)
wlam <- numeric(R)
wphi <- numeric(R)
wlamA1 <- numeric(R)
wphiA1 <- numeric(R)
wlamA2 <- numeric(R)
wphiA2 <- numeric(R)
wlamW <- numeric(R)
eps <- 10^(-6)

# start the simulation
for (r in 1:R)
  {
    cat("run",r,"\n")

    # generate mu
    mu <- rep(10,n)
    XX <- simfx1(n = n, p = p, tps = Tps)
    for (jj in 1:length(XX$funx))
      {
        mu <- mu + XX$funx[[jj]]%*%btrue[[jj]]
      }
    
    # generate y
    y <- mu + rnorm(n,0,sqrt(sigma2))

    # cross-validation
    cvError <- cv.grplMFLM(k = 5, Y = y, X = XX$funcs, Tps = Tps, lambda = lambda,
    phi = phi, dfs = 30)
    cvError <- apply(cvError,c(2,3),sum)
    wphi[r] <- wp <- which.min(apply(cvError,1,min))
    wlam[r] <- wl <- which.min(cvError[wp,])

    # fitting
    grpl <- grplMFLM(Y = y, X = XX$funcs, Tps = Tps, lambda = lambda[wl], phi = phi[wp],
    dfs = 30)

    # fitted coefs
    intercepts[r] <- grpl$intercept
    for (j in 1:p)
      {
        EstCoefs[[j]][,r] <- grpl$Coef[[j]]
      }
      
    # refund fit
    smoothFit <- pfr_old(Y = y, funcs = XX$funcs, kb = min(30,ceiling(n/p) - 1))
    interceptsS[r] <- smoothFit$beta.covariates

    # weights for adaptive solution/version 1
    wf1 <- wf2 <- numeric(p)
    for (j in 1:p)
      {
        bf <- EstCoefsSmooth[[j]][,r] <- smoothFit$BetaHat[[j]]
        wf1[j] <- 1/sqrt(sum(bf^2))
        tmax <- length(Tps[[j]])
        wf2[j] <- 1/(sqrt(sum((bf[-c(tmax-1,tmax)] - 2*bf[-c(1,tmax)] + bf[-c(1,2)])^2)) + eps)
        # second derivative of beta: 1/norm(second derivative)
      }
    wf1 <- p*wf1/sum(wf1)
    wf2 <- p*wf2/sum(wf2)

    # cross-validation
    cvError <- cv.grplMFLM(k = 5, Y = y, X = XX$funcs, Tps = Tps, lambda = lambda,
    phi = phi, dfs = 30, adapt1 = wf1, adapt2 = wf2)
    cvError <- apply(cvError,c(2,3),sum)
    wphiA1[r] <- wp <- which.min(apply(cvError,1,min))
    wlamA1[r] <- wl <- which.min(cvError[wp,])

    # fitting
    grpl <- grplMFLM(Y = y, X = XX$funcs, Tps = Tps, lambda = lambda[wl], phi = phi[wp],
    dfs = 30, adapt1 = wf1, adapt2 = wf2)

    # fitted coefs
    interceptsA1[r] <- grpl$intercept
    for (j in 1:p)
      {
        EstCoefsAdapt1[[j]][,r] <- grpl$Coef[[j]]
      }
    
    # weights for adaptive solution/version 2
    wf2 <- rep(1,length(wf2))

    # cross-validation
    cvError <- cv.grplMFLM(k = 5, Y = y, X = XX$funcs, Tps = Tps, lambda = lambda,
    phi = phi, dfs = 30, adapt1 = wf1, adapt2 = wf2)
    cvError <- apply(cvError,c(2,3),sum)
    wphiA2[r] <- wp <- which.min(apply(cvError,1,min))
    wlamA2[r] <- wl <- which.min(cvError[wp,])

    # fitting
    grpl <- grplMFLM(Y = y, X = XX$funcs, Tps = Tps, lambda = lambda[wl], phi = phi[wp],
    dfs = 30, adapt1 = wf1, adapt2 = wf2)

    # fitted coefs
    interceptsA2[r] <- grpl$intercept
    for (j in 1:p)
      {
        EstCoefsAdapt2[[j]][,r] <- grpl$Coef[[j]]
      }
    
    # cross-validation without smoothing
    cvError <- cv.grplMFLM(k = 5, Y = y, X = XX$funcs, Tps = Tps, lambda = lambda,
                           phi = 0, dfs = 30)
    cvError <- apply(cvError,c(2,3),sum)
    wlamW[r] <- wl <- which.min(cvError[1,])
    
    # fitting
    grpl <- grplMFLM(Y = y, X = XX$funcs, Tps = Tps, lambda = lambda[wl], phi = 0,
                     dfs = 30)
    
    # fitted coefs
    interceptsW[r] <- grpl$intercept
    for (j in 1:p)
    {
      EstCoefsWiggly[[j]][,r] <- grpl$Coef[[j]]
    }
  }


#save.image("sim2.RData")

# saved workspace

#postscript("PaperFGRPL/sim2stdrd.ps", width=12, height=8)
par(mfrow=c(2,5))
for (j in 1:10)
  {
    plot(Tps[[j]],btrue[[j]], type="l", col=1, ylim = c(-0.1,0.5),#range(EstCoefs[[j]]),
    lty=2,lwd=2, xlab="t", ylab=expression(beta(t)))
    #for (r in 1:R)
    #  lines(Tps[[j]],EstCoefs[[j]][,r])
    lines(Tps[[j]],apply(EstCoefs[[j]],1,mean),col=1,lwd=2)
    lines(Tps[[j]],apply(EstCoefs[[j]],1,mean) + apply(EstCoefs[[j]],1,sd),
    col=1,lty=3) #confidence interval?
    lines(Tps[[j]],apply(EstCoefs[[j]],1,mean) - apply(EstCoefs[[j]],1,sd),
    col=1,lty=3)
    title(paste("curve",j))
  }
#dev.off()

#postscript("PaperFGRPL/sim2adapt1.ps", width=12, height=8)
# par(mfrow=c(2,5))
# for (j in 1:10)
#   {
#     plot(Tps[[j]],btrue[[j]], type="l", col=1, ylim = c(-0.1,0.5),#range(EstCoefsAdapt[[j]]),
#     lty=2,lwd=2, xlab="t", ylab=expression(beta(t)))
#     #for (r in 1:R)
#     #  lines(Tps[[j]],EstCoefs[[j]][,r])
#     lines(Tps[[j]],apply(EstCoefsAdapt1[[j]],1,mean),col=1,lwd=2)
#     lines(Tps[[j]],apply(EstCoefsAdapt1[[j]],1,mean) + apply(EstCoefsAdapt1[[j]],1,sd),
#     col=1,lty=3)
#     lines(Tps[[j]],apply(EstCoefsAdapt1[[j]],1,mean) - apply(EstCoefsAdapt1[[j]],1,sd),
#     col=1,lty=3)
#     title(paste("curve",j))
#   }
#dev.off()

#postscript("PaperFGRPL/sim2adapt2.ps", width=12, height=8)
# par(mfrow=c(2,5))
# for (j in 1:10)
#   {
#     plot(Tps[[j]],btrue[[j]], type="l", col=1, ylim = c(-0.1,0.5),#range(EstCoefsAdapt[[j]]),
#     lty=2,lwd=2, xlab="t", ylab=expression(beta(t)))
#     #for (r in 1:R)
#     #  lines(Tps[[j]],EstCoefs[[j]][,r])
#     lines(Tps[[j]],apply(EstCoefsAdapt2[[j]],1,mean),col=1,lwd=2)
#     lines(Tps[[j]],apply(EstCoefsAdapt2[[j]],1,mean) + apply(EstCoefsAdapt2[[j]],1,sd),
#     col=1,lty=3)
#     lines(Tps[[j]],apply(EstCoefsAdapt2[[j]],1,mean) - apply(EstCoefsAdapt2[[j]],1,sd),
#     col=1,lty=3)
#     title(paste("curve",j))
#   }
#dev.off()

#postscript("PaperFGRPL/sim2smooth.ps", width=12, height=8)
# par(mfrow=c(2,5))
# for (j in 1:10)
#   {
#     plot(Tps[[j]],btrue[[j]], type="l", col=1, ylim = c(-0.1,0.5),#range(EstCoefsSmooth[[j]]),
#     lty=2,lwd=2, xlab="t", ylab=expression(beta(t)))
#     #for (r in 1:R)
#     #  lines(Tps[[j]],EstCoefs[[j]][,r])
#     lines(Tps[[j]],apply(EstCoefsSmooth[[j]],1,mean),col=1,lwd=2)
#     lines(Tps[[j]],apply(EstCoefsSmooth[[j]],1,mean) + apply(EstCoefsSmooth[[j]],1,sd),
#     col=1,lty=3)
#     lines(Tps[[j]],apply(EstCoefsSmooth[[j]],1,mean) - apply(EstCoefsSmooth[[j]],1,sd),
#     col=1,lty=3)
#     title(paste("curve",j))
#   }
#dev.off()

#postscript("PaperFGRPL/sim2simple.ps", width=12, height=8)
# par(mfrow=c(2,5))
# for (j in 1:10)
# {
#   plot(Tps[[j]],btrue[[j]], type="l", col=1, ylim = c(-0.1,0.5),#range(EstCoefsSmooth[[j]]),
#        lty=2,lwd=2, xlab="t", ylab=expression(beta(t)))
#   #for (r in 1:R)
#   #  lines(Tps[[j]],EstCoefs[[j]][,r])
#   lines(Tps[[j]],apply(EstCoefsWiggly[[j]],1,mean),col=1,lwd=2)
#   lines(Tps[[j]],apply(EstCoefsWiggly[[j]],1,mean) + apply(EstCoefsWiggly[[j]],1,sd),
#         col=1,lty=3)
#   lines(Tps[[j]],apply(EstCoefsWiggly[[j]],1,mean) - apply(EstCoefsWiggly[[j]],1,sd),
#         col=1,lty=3)
#   title(paste("curve",j))
# }
#dev.off()

# number of runs considered
R <- 50
#R <- 20

# Selection
selection <- matrix(0,R,length(EstCoefs))
selectionA1 <- matrix(0,R,length(EstCoefsAdapt1))
selectionA2 <- matrix(0,R,length(EstCoefsAdapt2))
selectionW <- matrix(0,R,length(EstCoefsWiggly))
for (j in 1:length(EstCoefs))
  {
    selection[,j] <- as.numeric(apply(EstCoefs[[j]]^2,2,sum) > 0)[1:R]
    selectionA1[,j] <- as.numeric(apply(EstCoefsAdapt1[[j]]^2,2,sum) > 0)[1:R]
    selectionA2[,j] <- as.numeric(apply(EstCoefsAdapt2[[j]]^2,2,sum) > 0)[1:R]
    selectionW[,j] <- as.numeric(apply(EstCoefsWiggly[[j]]^2,2,sum) > 0)[1:R]
  }
size <- apply(selection,1,sum)
sizeA1 <- apply(selectionA1,1,sum)
sizeA2 <- apply(selectionA2,1,sum)
sizeW <- apply(selectionW,1,sum)
avSize <- c(mean(size),mean(sizeA1),mean(sizeA2),mean(sizeW))
avSize

Select <- rbind(apply(selection,2,mean),
                apply(selectionA1,2,mean),
                apply(selectionA2,2,mean),
                apply(selectionW,2,mean))
Select <- round(Select, digits=3)
Select <- cbind(Select,avSize)
row.names(Select) <- c("stdrd","adapt1","adapt2","simple")
#library(xtable)
#xSc <- xtable(Select, digits=2)
#print(xSc, file="PaperFGRPL/sim2select.tex")

# xnames <- c("stdrd","adapt1","adapt2","simple")
# par(mfrow=c(2,5))
# for (j in 1:10)
#   {
#     barplot(c(apply(selection,2,mean)[j],
#               apply(selectionA1,2,mean)[j],
#               apply(selectionA2,2,mean)[j],
#               apply(selectionW,2,mean)[j]), names.arg=xnames, ylim=c(0,1))
#   }
  

mse <- matrix(0,R,length(EstCoefs))
mseA1 <- matrix(0,R,length(EstCoefsAdapt1))
mseA2 <- matrix(0,R,length(EstCoefsAdapt2))
mseW <- matrix(0,R,length(EstCoefsWiggly))
mseS <- matrix(0,R,length(EstCoefsSmooth))
for (j in 1:length(EstCoefs))
  {
    mse[,j] <- apply((EstCoefs[[j]] - btrue[[j]])^2,2,sum)
    mseA1[,j] <- apply((EstCoefsAdapt1[[j]] - btrue[[j]])^2,2,sum)
    mseA2[,j] <- apply((EstCoefsAdapt2[[j]] - btrue[[j]])^2,2,sum)
    mseW[,j] <- apply((EstCoefsWiggly[[j]] - btrue[[j]])^2,2,sum)
    mseS[,j] <- apply((EstCoefsSmooth[[j]] - btrue[[j]])^2,2,sum)    
  }
apply(mse,2,mean)
apply(mseA1,2,mean)
apply(mseA2,2,mean)
apply(mseW,2,mean)
apply(mseS,2,mean)


# xnames <- c("pfr","simple","adapt2","adapt1","stdrd")
# ylims <- cbind(rep(0,6),c(7,3,1,9,5,2))
# #postscript("PaperFGRPL/sim2mse.ps", width=8, height=2.5, horizontal=F)
# par(mfrow=c(1,6),bty="n", mar=c(1, 4, 4, 0))
# boxplot(cbind(mse[,1],mseA1[,1],mseA2[,1],mseW[,1],mseS[,1]), col=c(2,3,4,5,6),
# names=NULL, horizontal=F, las=0, ylab="SE (dense design, normal y)",
# xaxt="n", ylim=ylims[1,])
# #title(paste("curve 1"))
# for (j in 2:5)
#   {
#     boxplot(cbind(mse[,j],mseA1[,j],mseA2[,j],mseW[,j],mseS[,j]), col=c(2,3,4,5,6),
#     names=NULL, horizontal=F, las=0, ylab=" ", xaxt="n", ylim=ylims[j,])
#     #title(paste("curve",j))
#   }
# boxplot(cbind(c(mse[,6:10]),c(mseA1[,6:10]),c(mseA2[,6:10]),
# c(mseW[,6:10]),c(mseS[,6:10])), col=c(2,3,4,5,6), names=NULL,
# horizontal=F, las=0, xaxt="n", ylab=" ", ylim=ylims[6,])
# #title("curve 6-10")
# #dev.off()
# 
# 
# #postscript("PaperFGRPL/sim2curves.ps", width=8, height=0.9, horizontal=F)
# par(mfrow=c(1,6),bty="n", mar=c(1, 4, 2, 0))
# plot(btrue[[1]],type="l", xaxt="n", xlab=" ", ylab=" ", ylim=c(0,0.49), lwd=2)
# title("curve 1")
# for (j in 2:5)
#   {
#     plot(btrue[[j]],type="l", axes=F, xlab=" ", ylab=" ", ylim=c(0,0.49), lwd=2)
#     title(paste("curve",j))
#   }
# plot(btrue[[6]],type="l", axes=F, xlab=" ", ylab=" ", ylim=c(0,0.49), lwd=2)
# title("curve 6-10")
# #dev.off()
# 
# 
# 


# test data set
set.seed(1111)
m <- 5000
mu <- rep(10,m)
XX <- simfx1(n = m, p = p, tps = Tps)
for (jj in 1:length(XX$funx))
  {
    mu <- mu + XX$funx[[jj]]%*%btrue[[jj]]
  }
y <- mu + rnorm(m,0,sqrt(sigma2))


# prediction
msep <- numeric(R)
msepA1 <- numeric(R)
msepA2 <- numeric(R)
msepS <- numeric(R)
msepW <- numeric(R)
for (r in 1:R)
  {
    muEst <- rep(intercepts[r],m)
    muEstA1 <- rep(interceptsA1[r],m)
    muEstA2 <- rep(interceptsA2[r],m)
    muEstS <- rep(interceptsS[r],m)
    muEstW <- rep(interceptsW[r],m)
    for (jj in 1:length(XX$funx))
      {
        muEst <- muEst + XX$funcs[[jj]]%*%EstCoefs[[jj]][,r]
        muEstA1 <- muEstA1 + XX$funcs[[jj]]%*%EstCoefsAdapt1[[jj]][,r]
        muEstA2 <- muEstA2 + XX$funcs[[jj]]%*%EstCoefsAdapt2[[jj]][,r]
        muEstS <- muEstS + XX$funcs[[jj]]%*%EstCoefsSmooth[[jj]][,r]
        muEstW <- muEstW + XX$funcs[[jj]]%*%EstCoefsWiggly[[jj]][,r]
      }
    msep[r] <- mean((y - muEst)^2)
    msepA1[r] <- mean((y - muEstA1)^2)
    msepA2[r] <- mean((y - muEstA2)^2)
    msepS[r] <- mean((y - muEstS)^2)
    msepW[r] <- mean((y - muEstW)^2)
  }


# xnames <- c("pfr","simple","adapt2","adapt1","stdrd")
# #postscript("PaperFGRPL/sim2msep.ps", width=3, height=4, horizontal=F)
# par(mfrow=c(1,1), bty="n", mar=c(1, 4, 4, 0))
# boxplot(cbind(msep,msepA1,msepA2,msepW,msepS), names=NULL, col=c(2,3,4,5,6),
# xaxt="n", horizontal=F, las=0, ylab="mean squared error of prediction",
# ylim=c(10.4,12.8))
# title("dense design, normal y")
# #dev.off()
save.image("sim2_50runs.RData")