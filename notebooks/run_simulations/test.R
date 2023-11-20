library(here)

simulation_name = "3_predictors"
simulation_settings_file = "default"
# Source the generic simulator script
source(here("src", "R", "generic_simulator", "simulate_test.R"))
source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))


# Required inputs before running the simulation!!!
inputs  <- load_simulation_settings(simulation_name, simulation_settings_file)
time_domains_eval <- lapply(inputs$time_domains, function(domain) {
           seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
       })
       inputs$time_domains <- time_domains_eval
inputs$noise_snr = c(25,25)

betas = create_beta_curves(beta_funcs, inputs$time_domains)

X <- simulate_functional_features(mu_funcs, cov_funcs, 1, inputs$time_domains)

Y = compute_Y_values(X, betas, observations = 1, predictors, time_domains = inputs$time_domains, intercept=inputs$intercept)