# config.R for setting utility paths relative to the project root directory
library(here)

# Source utility files using the relative paths
source(here("src", "R", "data_simulator",  "utils", "basis_utilities.R"))
source(here("src", "R", "data_simulator",  "utils", "covariance_utilities.R"))
source(here("src", "R", "data_simulator",  "utils", "simulation_utilities.R"))
source(here("src", "R", "data_simulator",  "utils", "model_utilities.R"))
source(here("src", "R", "data_simulator",  "utils", "plot_utilities.R"))
