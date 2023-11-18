library(here)

simulation_name = "paper"
simulation_settings_file = "default"
source(here("simulations","load_and_run.R"))

print(inputs$observations)
time_domains_eval <- lapply(inputs$time_domains, function(domain) {
    seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
})
inputs$time_domains <- time_domains_eval
outputs <- run_simulation(inputs)