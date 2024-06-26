# ------------------------------------------------------------------------------
# -------- Functions for Generating Functional Data Sets -----------------------
# ------------------------------------------------------------------------------



# -------- generate functional predictors directly -----------------------------
# -------- by assuming a specific form             -----------------------------

# similar to Tutz & Gertheiss (2010)
simfx1 <- function(n, p, tps, varx=rep(0,p), bx=5, mx=2*pi)
  {
    fx <- list()
    fobs <- list()
    for (j in 1:p)
      {
        tmax <- max(tps[[j]])
        fx[[j]] <- matrix(0, n, length(tps[[j]]))
        fobs[[j]] <- matrix(0, n, length(tps[[j]]))
        for (i in 1:n)
          {
            bij <- runif(5, 0, bx)
            mij <- runif(5, 0, mx)
            tfx <- function(tp)
              {
                (sum(bij*sin(tp*(5-bij)*(2*pi/tmax)) - mij) + 15)/100
              }
            fx[[j]][i,] <- sapply(tps[[j]],tfx)
          }
      }
    fx <- lapply(fx, scale)
    for (j in 1:p)
    {
      fx[[j]] <- fx[[j]]/10
      for (i in 1:n)
      {
        fobs[[j]][i,] <- fx[[j]][i,] + rnorm(length(tps[[j]]), 0, sqrt(varx[j]))        
      }
    }
    return(list("funx"=fx, "funcs"=fobs))
  }
