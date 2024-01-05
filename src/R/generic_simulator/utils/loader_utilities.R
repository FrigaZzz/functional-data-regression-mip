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
        variant_path <- here("simulations", "settings", setting_name, "variants", paste0(variant_name, ".R"))
        if (file.exists(variant_path)) {
            source(variant_path)
        } else {
            warning("Variant settings file does not exist: ", variant_path)
        }
    }

    # check if true_predictors exists
    if (!exists("simulation_type")) {
        simulation_type <- "non paper"
    }
    # Usage:
    params <- list(
        predictors = predictors,
        observations = observations,
        measurements = measurements,
        basis_functions = basis_functions,
        intercept = intercept,
        norder = norder,
        mu_funcs = mu_funcs,
        beta_funcs = beta_funcs,
        time_domains = time_domains,
        cov_funcs = cov_funcs,
        seed = seed,
        noise_snr = noise_snr,
        true_predictors = true_predictors,
        simulation_type = simulation_type
    )
    # Return the environment where the settings are stored
    return(params)
}


