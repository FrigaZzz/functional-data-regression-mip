using RCall


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
function load_simulation_data(simulation_name, simulation_settings_file, project_root; observations=nothing, measurements=nothing, basis_functions=nothing, error_sd=nothing, seed=nothing)
    # clear R environment
    # load new simulation data
    runner_file_path = joinpath(project_root, "simulations", "load_and_run.R")

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
    if error_sd != nothing
        overrides["error_sd"] = error_sd
    end
    if seed != nothing
        overrides["seed"] = seed
    end

    @rput simulation_name
    @rput simulation_settings_file
    @rput runner_file_path
    @rput overrides
    R"""
        source(runner_file_path)
        for (param_name in names(overrides)) {
            inputs[[param_name]] <- overrides[[param_name]]
        }
        print(inputs$observations)
        outputs <- run_simulation(inputs)
    """
    full_output = rcopy(R"(outputs)")
    full_input = rcopy(R"(inputs)")
    output = merge(full_input, full_output)
    return output
end