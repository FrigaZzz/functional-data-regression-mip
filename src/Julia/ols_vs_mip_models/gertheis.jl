using JuMP
using Gurobi
using Statistics
using RCall



function mip_functional_regression(
                                    Y, X, 
                                    UpperBounds, LowerBounds; 
                                    intercept = false, 
                                    group_limit=Inf, lambda = 0,
                                    initial_beta = nothing, initial_group = nothing, initial_alpha = nothing)

    @rput Y X  # Pass Julia variables to R

    R"""

    library(here)
    source(here("src", "Julia",  "ols_vs_mip_models",   "gerth.R"))

    
    obseravation <- dim(X)[1]    
    predictors <- dim(X)[2]
    measurements <- dim(X)[3]

    fobs <- list()
    for (j in 1:predictors) {
      fobs[[j]] <- matrix(0, obseravation, measurements)
      for (i in 1:obseravation) {
        fobs[[j]][i,] <- X[i,j,]
      }
    }

    tps <- list()
    # for any time_domain, create a vector of time points from 1:measurements
    tps <- inputs$time_domains

    # fit using functional smooth group lasso
    lambda <- 10^seq(3,-3,by=-1)
    phi <- 10^seq(3, -3, by=-1)

    # cross-validation
    cvError <- grplMFLM(k = 2, Y = Y, X = fobs, Tps = tps, lambda = lambda , phi = phi, dfs = inputs$basis_functions)
    cvError <- apply(cvError,c(2,3),sum)
    best_phi <- wp <- which.min(apply(cvError,1,min))
    best_lambda <- wl <- which.min(cvError[wp,])
    best_phi_ <- best_phi[1]
    best_lambda_ <- best_lambda
    grpl = grplFlinear(Y = Y, X = fobs, Tps = tps, lambda = lambda[best_lambda], phi = phi[best_phi_], dfs =inputs$basis_functions )
    intercept <- grpl$intercept 
    betas <- grpl$real_coef[-1, ]

    print(paste("lambda: ", lambda[best_lambda], " phi: ", phi[best_phi_]))

    matrix = do.call(rbind, lapply(grpl$Coef, t))
    ## convert the given betas into our format
    source(here("src","R", "generic_simulator", "config.R"))  
    result <- smooth_betas_generic(matrix, inputs$basis_functions, inputs$time_domains, outputs$basis_objs)
    beta_out = result$beta_matrix

    #beta_out <- (to_matrix_form(betas,predictors = predictors, basis_functions = inputs$basis_functions))
    """

    @rget beta_out intercept  # Get the R variable back into Julia

    beta_star = beta_out
    group_star = ones(Int, predictors)  # Create a vector of 1s with length predictors
    a_star = intercept
    return beta_star, a_star, group_star
end