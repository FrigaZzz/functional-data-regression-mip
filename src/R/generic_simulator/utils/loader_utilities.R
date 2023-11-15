library(here)

#' Load simulation settings with optional variant override
#'
#' This function loads simulation settings from a default file and, optionally, from a variant file.
#' The default file must exist, otherwise an error is thrown.
#' If a variant file is specified but does not exist, a warning is issued.
#' The function returns the environment where the settings are stored.
#'
#' @param setting_name A character string specifying the name of the simulation settings.
#' @param variant_name A character string specifying the name of the variant settings, if any.
#'
#' @return The environment where the settings are stored.
#'
#' @examples
#' # Assume that 'default.R' defines variables a, b, c, and 'variant1.R' overrides variable b
#' settings_env <- load_simulation_settings("sim1", "var1")
#' a <- settings_env$a
#' b <- settings_env$b
#' c <- settings_env$c
load_simulation_settings <- function(setting_name, variant_name = NULL) {
    default_path <- here("simulations", "settings", setting_name, "default.R")
    variant_path <- NULL
    
    # Load the default settings
    if (file.exists(default_path)) {
        source(default_path)
    } else {
        stop("Default settings file does not exist: ", default_path)
    }
    
    # If a variant is specified, construct its path and source it
    if (!is.null(variant_name) && variant_name != "default") {
        variant_path <- here( "simulations", "settings", setting_name, "variants", paste0(variant_name, ".R"))
        if (file.exists(variant_path)) {
            source(variant_path)
        } else {
            warning("Variant settings file does not exist: ", variant_path)
        }
    }
    # check if true_predictors exists
    if (!exists("true_predictors")) {
        true_predictors <- non_zero_betas(beta_funcs)
    }
    # Usage:
    params <- list(
        predictors = predictors,
        measurements = measurements,
        observations = observations,
        basis_functions = basis_functions,
        intercept = intercept,
        norder = norder,
        error_sd = error_sd,
        seed = seed,
        noise_sd = noise_sd,
        mu_funcs = mu_funcs,
        cov_funcs = cov_funcs,
        beta_funcs = beta_funcs,
        time_domains = time_domains,
        true_predictors = true_predictors
    )
    # Return the environment where the settings are stored
    return(params)
}


#' Count the number of non-zero functions in a list of functions
#'
#' This function takes a list of functions and counts the number of non-zero functions.
#' 
#' @param func_list A list of functions to be counted.
#' 
#' @return An integer indicating the number of non-zero functions in the list.
#'
#' @examples
#' count_non_zero_functions(list(function(x) x^2, function(x) 0*x))
#' # returns 1
#'
non_zero_betas<- function(func_list) {
    zero_funcs <- sapply(func_list, function(f) {
        body_text <- deparse(body(f))
        is.zero <- grepl("0 \\*", body_text) 
        # return 0 if true, 1 if false
        return(1 - is.zero)
    })
    return(zero_funcs)
}
