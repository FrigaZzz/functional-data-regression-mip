c_val <- 0.4
true_predictors <- c(1,0,0,0,0,0)


beta_funcs = list(
  function(t) {
    sin(t)
  },
  function(t) {
    rep(0, length(t))  
  },
  function(t) {
    rep(0, length(t))  
  },
  function(t) {
    rep(0, length(t))  
  },
  function(t) {
    rep(0, length(t))  
  },
  function(t) {
    rep(0, length(t))  # This function always returns 0 regardless of the input t
  }
)

