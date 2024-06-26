using RCall
R"""
library(here)
"""

"""
    load_simulation_data(simulation_name, simulation_settings_file, project_root; observations=nothing, measurements=nothing, basis_functions=nothing, error_sd=nothing, seed=nothing)

Function to load R simulation data.

# Arguments
- `simulation_name::String`: Name of the simulation.
- `simulation_settings_file::String`: Path to the simulation settings file.
- `project_root::String`: Path to the project root directory.
- `observations::Union{Nothing, Array{Float64, 2}}=nothing`: Observations matrix.
- `measurements::Union{Nothing, Array{Float64, 2}}=nothing`: Measurements matrix.
- `basis_functions::Union{Nothing, Array{Float64, 2}}=nothing`: Basis functions matrix.
- `error_sd::Union{Nothing, Float64}=nothing`: Standard deviation of the error.
- `seed::Union{Nothing, Int64}=nothing`: Seed for the random number generator.

# Returns
- `output::Dict`: A dictionary containing the simulation inputs and outputs.
"""
# Function to load R simulation data
function load_simulation_data(simulation_name, simulation_settings_file, project_root; observations=nothing, measurements=nothing, basis_functions=nothing, noise_snr=nothing, seed=nothing)
    # clear R environment
    # load new simulation data
    overrides = Dict()
    if observations != nothing
        overrides["observations"] = observations
    end
    if measurements != nothing
        overrides["measurements"] = measurements
    end
    if basis_functions != nothing
        overrides["basis_functions"] = basis_functions
    end
    if noise_snr != nothing
        overrides["noise_snr"] = noise_snr
    end
    if seed != nothing
        overrides["seed"] = seed
    end

    @rput simulation_name
    @rput simulation_settings_file
    @rput overrides
    R"""
    # Source the generic simulator script
        source(here("src", "R", "generic_simulator", "simulate_main.R"))
        source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))

        inputs  <- load_simulation_settings(simulation_name, simulation_settings_file)

        for (param_name in names(overrides)) {
            inputs[[param_name]] <- overrides[[param_name]]
        }
        time_domains_eval <- lapply(inputs$time_domains, function(domain) {
            seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
        })
        inputs$time_domains <- time_domains_eval

      

        outputs <- result <- do.call(generate_data, inputs[c("predictors", "observations", "measurements", "basis_functions", "intercept", "norder", "mu_funcs", "beta_funcs","time_domains", "cov_funcs", "seed","noise_snr","simulation_type")] )

    """
    full_output = rcopy(R"(outputs)")
    full_input = rcopy(R"(inputs)")
    output = merge(full_input, full_output)
    return output
end

# Function to load R simulation data
function load_simulation_reiman(simulation_name, simulation_settings_file, project_root; observations=nothing, measurements=nothing, basis_functions=nothing, noise_snr=nothing, seed=nothing)
    overrides = Dict()
    if observations != nothing
        overrides["observations"] = observations
    end
    if measurements != nothing
        overrides["measurements"] = measurements
    end
    if basis_functions != nothing
        overrides["basis_functions"] = basis_functions
    end
    if noise_snr != nothing
        overrides["noise_snr"] = noise_snr
    end
    if seed != nothing
        overrides["seed"] = seed
    end

    @rput simulation_name
    @rput simulation_settings_file
    @rput overrides
    R"""
        # Source the generic simulator script
        source(here("src", "R", "generic_simulator", "simulate_main_reiman.R"))
        source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))
        set.seed(overrides$seed)

        # Required inputs before running the simulation!!!
        inputs  <- load_simulation_settings(simulation_name,        simulation_settings_file)       
        for (param_name in names(overrides)) {
            inputs[[param_name]] <- overrides[[param_name]]
        }

        time_domains_eval <- lapply(inputs$time_domains, function(domain) {
            seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
        })
        inputs$time_domains <- time_domains_eval
        
        

        outputs <-  do.call(generate_data, inputs[c("predictors", "observations", "measurements", "basis_functions", "intercept", "norder", "mu_funcs", "beta_funcs","time_domains", "cov_funcs", "seed","noise_snr","simulation_type")] )
    """
    full_output = rcopy(R"(outputs)")
    full_input = rcopy(R"(inputs)")
    output = merge(full_input, full_output)
    return output
end


# Function to load R simulation data
function load_simulation_robust(simulation_name, simulation_settings_file, project_root; observations=nothing, measurements=nothing, basis_functions=nothing, noise_snr=nothing, seed=nothing)
    # clear R environment
    # load new simulation data
    overrides = Dict()
    if observations != nothing
        overrides["observations"] = observations
    end
    if measurements != nothing
        overrides["measurements"] = measurements
    end
    if basis_functions != nothing
        overrides["basis_functions"] = basis_functions
    end
    if noise_snr != nothing
        overrides["noise_snr"] = noise_snr
    end
    if seed != nothing
        overrides["seed"] = seed
    end

    @rput simulation_name
    @rput simulation_settings_file
    @rput overrides
    R"""
    # Source the generic simulator script
        source(here("src", "R", "generic_simulator", "simulate_main_robust.R"))
        source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))

        inputs  <- load_simulation_settings(simulation_name, simulation_settings_file)

        for (param_name in names(overrides)) {
            inputs[[param_name]] <- overrides[[param_name]]
        }
        time_domains_eval <- lapply(inputs$time_domains, function(domain) {
            seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
        })
        inputs$time_domains <- time_domains_eval

        outputs <-   result <- do.call(generate_data, inputs[c("predictors", "observations", "measurements", "basis_functions", "intercept", "norder", "mu_funcs", "beta_funcs","time_domains", "cov_funcs", "seed","noise_snr","simulation_type")] )
    """
    full_output = rcopy(R"(outputs)")
    full_input = rcopy(R"(inputs)")
    output = merge(full_input, full_output)
    return output
end


# Function to load R simulation data
function load_simulation_gertheiss(simulation_name, simulation_settings_file, project_root; observations=nothing, measurements=nothing, basis_functions=nothing, noise_snr=nothing, seed=nothing, phi=nothing)
    # clear R environment
    # load new simulation data

    overrides = Dict()
    if observations != nothing
        overrides["observations"] = observations
    end
    if measurements != nothing
        overrides["measurements"] = measurements
    end
    if basis_functions != nothing
        overrides["basis_functions"] = basis_functions
    end
    if noise_snr != nothing
        overrides["noise_snr"] = noise_snr
    end
    if seed != nothing
        overrides["seed"] = seed
    end

    if phi != nothing
        overrides["phi"] = phi
    end
    
    @rput simulation_name
    @rput simulation_settings_file
    @rput overrides
    R"""
        # Source the generic simulator script
        source(here("src", "R", "generic_simulator", "simulate_main_gertheiss.R"))
        source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))

        inputs  <- load_simulation_settings(simulation_name, simulation_settings_file)

        for (param_name in names(overrides)) {
            inputs[[param_name]] <- overrides[[param_name]]
        }
        time_domains_eval <- lapply(inputs$time_domains, function(domain) {
            seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
        })
        inputs$time_domains <- time_domains_eval
        print(inputs)
        # Required inputs before running the simulation!!!

        outputs <- do.call(generate_data, inputs[c("predictors", "observations", "measurements", "basis_functions", "intercept", "norder", "mu_funcs", "beta_funcs","time_domains", "cov_funcs", "seed","noise_snr","simulation_type","phi")] )

    """
    full_output = rcopy(R"(outputs)")
    full_input = rcopy(R"(inputs)")
    output = merge(full_input, full_output)
    return output
end