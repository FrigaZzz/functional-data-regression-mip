# config.R for setting utility paths relative to the project root directory
library(here)

# Source utility files using the relative paths
source(here("src", "R",  "generic_simulator",   "utils" , "covariance_utilities.R"))
source(here("src", "R",  "generic_simulator",   "utils" , "plot_utilities.R"))
source(here("src", "R",   "generic_simulator",  "utils" , "model_utilities.R")) # sets the utility 
source(here("src", "R",   "generic_simulator",  "utils" , "simulation_utilities.R")) # sets the utility paths
source(here("src", "R",   "generic_simulator",  "utils" , "basis_utilities.R")) # sets the utility paths