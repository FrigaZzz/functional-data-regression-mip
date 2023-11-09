
using Pkg
project_root = dirname(@__DIR__)

Pkg.activate(project_root)      # Activate the environment in the current project directory
Pkg.instantiate()      # Install all the dependencies as specified in the Project.toml and Manifest.toml
