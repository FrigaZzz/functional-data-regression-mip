# parameters_3_predictors.R
# Parameter definitions for a simulation with 3 predictors.

true_predictors <- c(1,0,0)

beta_funcs <- list(
  function(t) cos(t),
  function(t) 0 * t,
  function(t) 0 * t
)

