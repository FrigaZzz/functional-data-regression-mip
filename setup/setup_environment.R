# Initialize renv if it's not already present and restore the project environment
if (!require(renv)) {
  install.packages("renv")
}

# Install and load the here package if it's not already present
if (!require(here)) {
  install.packages("here")
  library(here)
}

# Change the working directory to the parent directory
setwd(here("."))

renv::init(bare = TRUE, restart = TRUE)
renv::restore()          # Restore packages from the renv.lock file