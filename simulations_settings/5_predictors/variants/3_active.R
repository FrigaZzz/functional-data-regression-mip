# parameters_5_predictors.R
# Parameter definitions for a simulation with 5 predictors.

true_predictors <- c(1,1,0,0,1)

beta_funcs <- list(
  function(t) 4 * exp(-abs(t - 1)),
  function(t) 2 * (1 - abs(t - 1)),
  function(t) 0 * t,
  function(t) 0 * t,
  function(t) 2 * t
)

