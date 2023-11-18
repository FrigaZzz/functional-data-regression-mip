
# Parameter definitions for a custom simulation
# vector of affected betas by c value
c_val <- 0.4
true_predictors <- c(1,1,1,1,1,0)


beta_funcs = list(
  function(t) {
    sin(t)
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    -c_val * t^2
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    c_val * sin(pi * t)
  },
  function(t) {
    0 *t # This function always returns 0 regardless of the input t
  }
)

