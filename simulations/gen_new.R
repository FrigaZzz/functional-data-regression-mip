

# Parameter definitions for a custom simulation
# vector of affected betas by c value
c_val <- 0
true_predictors <- c(1, 1, 0, 1, 0, 0)

predictors <- 6
measurements <- 1000
observations <- 250
basis_functions <- 6
intercept <- 0
norder <- 4
noise_snr = c(100,100)
seed <- 1

simulation_type = "paper"
beta_funcs <- list(
  function(t) {
    sin(t)
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    if (c_val == 0) {
      rep(0, length(t))
    } else {
      -c_val * t^2
    }
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    if (c_val == 0) {
      rep(0, length(t))
    } else {
      c_val * sin(pi * t)
    }
  },
  function(t) {
    rep(0, length(t)) # This function always returns 0 regardless of the input t
  }
)


mu_funcs <- list(
  function(t, args) {
    cos(2 * pi * (t - args$a1)) + args$a2
  },
  function(t, args) {
    args$b1 * sin(pi * t) + args$b2
  },
  function(t, args) {
    args$c1 * t^3 + args$c2 * t^2 + args$c3 * t
  },
  function(t, args) {
    sin(2 * (t - args$d1)) + args$d2 * t
  },
  function(t, args) {
    args$e1 * cos(2 * t) + args$e2 * t
  },
  function(t, args) {
    args$f1 * exp(-t / 3) + args$f2 * t + args$f3
  }
)

# time_domains <- list(
#   list(0, 1),
#   list(0, 1),
#   list(0, 1),
#   list(0, 1),
#   list(0, 1),
#   list(0, 1)
# )

time_domains <- list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1),
  list(0, pi / 3),
  list(-2, 1),
  list(-1, 1)
)
# cov_funcs is null only in the paper simulation
cov_funcs <- NULL


generate_coefficients <- function(observations, coef_specs) {
  coef_list <- list()
  for (predictor in names(coef_specs)) {
    predictor_specs <- coef_specs[[predictor]]
    predictor_coefs <- list()
    for (coef_name in names(predictor_specs)) {
      spec <- predictor_specs[[coef_name]]
      if (spec$type == "norm") {
        predictor_coefs[[coef_name]] <- rnorm(observations, mean = spec$mean, sd = spec$sd)
      } else if (spec$type == "unif") {
        predictor_coefs[[coef_name]] <- runif(observations, min = spec$min, max = spec$max)
      }
    }
    coef_list[[predictor]] <- predictor_coefs
  }
  return(coef_list)
}


coef_specs_original <- list(
  '1' = list(
    a1 = list(type = "norm", mean = -5, sd = 3),
    a2 = list(type = "norm", mean = 7, sd = 1.5)
  ),
  '2' = list(
    b1 = list(type = "unif", min = 4, max = 7),
    b2 = list(type = "norm", mean = 0, sd = 1)
  ),
  '3' = list(
    c1 = list(type = "norm", mean = -3, sd = 1.2),
    c2 = list(type = "norm", mean = 2, sd = 0.5),
    c3 = list(type = "norm", mean = -2, sd = 1),
    c4 = list(type = "norm", mean = 2, sd = 0.75) # only for the other paper simulation
  ),
  '4' = list(
    d1 = list(type = "norm", mean = -2, sd = 1),
    d2 = list(type = "norm", mean = 3, sd = 1.5)
  ),
  '5' = list(
    e1 = list(type = "unif", min = 4, max = 7),
    e2 = list(type = "norm", mean = 2, sd = 0.4)
  ),
  '6' = list(
    f1 = list(type = "norm", mean = 4, sd = 2),
    f2 = list(type = "norm", mean = -3, sd = 0.5),
    f3 = list(type = "norm", mean = 1, sd = 1)
  )
)



#' 
#' This function simulates functional features for the paper using the provided mean functions,  number of observations and time domains.
#' 
#' @param mu_funcs A list of mean functions.
#' @param observations The number of observations.
#' @param time_domains A list of time domains.
#' 
#' @return An array of simulated functional features.
#' 
simulate_true_predictors_Ut <- function(mu_funcs, observations, time_domains) {
  predictors = length(mu_funcs)
  measurements = length(time_domains[[1]])
  coef_list <- generate_coefficients(observations, coef_specs)


  X <- array(0, dim = c(observations, predictors, measurements))
  # Generate the functional data for each observation and each predictor
  for (i in 1:observations) {
    for (j in 1:predictors) {
      # Extract the coefficients for the j-th predictor of the i-th observation
      current_coefs <- lapply(coef_list[[j]], function(row) row[i])
      # Apply the mu_funcs function with the time domain for the j-th predictor
      # and the current coefficients for the i-th observation
      X[i, j, ] <- mu_funcs[[j]](time_domains[[j]], current_coefs)
    }
  }
  
  return(X)
}

generate_coefficients <- function(observations, coef_specs) {
  coef_list <- list()
  for (predictor in names(coef_specs)) {
    predictor_specs <- coef_specs[[predictor]]
    predictor_coefs <- list()
    for (coef_name in names(predictor_specs)) {
      spec <- predictor_specs[[coef_name]]
      if (spec$type == "norm") {
        predictor_coefs[[coef_name]] <- rnorm(observations, mean = spec$mean, sd = spec$sd)
      } else if (spec$type == "unif") {
        predictor_coefs[[coef_name]] <- runif(observations, min = spec$min, max = spec$max)
      }
    }
    coef_list[[predictor]] <- predictor_coefs
  }
  return(coef_list)
}


coef_specs_original <- list(
  '1' = list(
    a1 = list(type = "norm", mean = -5, sd = 3),
    a2 = list(type = "norm", mean = 7, sd = 1.5)
  ),
  '2' = list(
    b1 = list(type = "unif", min = 4, max = 7),
    b2 = list(type = "norm", mean = 0, sd = 1)
  ),
  '3' = list(
    c1 = list(type = "norm", mean = -3, sd = 1.2),
    c2 = list(type = "norm", mean = 2, sd = 0.5),
    c3 = list(type = "norm", mean = -2, sd = 1),
    c4 = list(type = "norm", mean = 2, sd = 0.75) # only for the other paper simulation
  ),
  '4' = list(
    d1 = list(type = "norm", mean = -2, sd = 1),
    d2 = list(type = "norm", mean = 3, sd = 1.5)
  ),
  '5' = list(
    e1 = list(type = "unif", min = 4, max = 7),
    e2 = list(type = "norm", mean = 2, sd = 0.4)
  ),
  '6' = list(
    f1 = list(type = "norm", mean = 4, sd = 2),
    f2 = list(type = "norm", mean = -3, sd = 0.5),
    f3 = list(type = "norm", mean = 1, sd = 1)
  )
)
# Assuming generate_coefficients and coef_specs_original are defined
coef_specs <- coef_specs_original 

# Run the simulation
simulated_data <- simulate_true_predictors_Ut(mu_funcs, observations, time_domains)

