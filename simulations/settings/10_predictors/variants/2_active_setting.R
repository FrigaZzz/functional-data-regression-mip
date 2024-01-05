# parameters_10_predictors.R
# This file defines parameters for a simulation with 10 predictors.


true_predictors <- c(1,1,0,0,0,0,0,0,0,0)

  beta_funcs = list(
    function(t) 4 * exp(-abs(t - 1)),
    function(t) 2 * (1 - abs(t - 1)),
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t, # This assumes that the function returns a single scalar zero
    function(t) 0 * t  # This assumes that the function returns a single scalar zero
  )
