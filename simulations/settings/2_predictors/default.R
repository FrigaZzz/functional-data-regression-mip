# parameters_3_predictors.R
# Parameter definitions for a simulation with 3 predictors.

beta_funcs <- list(
  function(t) (2 * t) + 1,
  function(t) rep(0, length(t)) 
)
true_predictors <- c(1,0)

predictors <- 2
measurements <- 50
observations <- 5000
basis_functions = 5
intercept = 0
norder = 4
error_sd = 0.01
noise_sd = 0.00
seed = 1

mu_funcs <- list(
   function(t,args) {
      cos( pi * (t - args$a1)) + args$a2
        
    },
    function(t,args) {

        args$b1 * sin(pi * t) + args$b2
    }
)

cov_funcs <- NULL

time_domains = list(
  list(0,  pi / 2),
  list(0, pi / 2)
)
