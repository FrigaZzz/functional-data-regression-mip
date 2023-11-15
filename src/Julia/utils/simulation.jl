using RCall


# Function to load R simulation data
function load_simulation_data(simulation_name, simulation_settings_file, project_root, observations, measurements, basis_functions, error_sd, seed)
    # clear R environment
    # load new simulation data

    simulation_file_path = joinpath(project_root, "simulations", "settings", simulation_name, simulation_settings_file * ".R")
    runner_file_path = joinpath(project_root, "simulations", "run.R")

    # Optionally, define overrides here and pass them to R
    overrides = Dict()
    # fill overrides with non null/nothing values from the function arguments
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
    @rput simulation_file_path
    @rput runner_file_path
    @rput overrides

    R"""
    source(simulation_file_path)

    # Apply overrides if provided
    if (exists('overrides')) {
        for (param_name in names(overrides)) {
            param_name <- overrides[[param_name]]
        }
    }
    
    source(runner_file_path)
    """
    # Load simulation data, setup directories, etc., now using the passed simulation_name
    # Setup variables from R environment
    # Sample input data
    predictors = rcopy(R"(params$predictors)")
    true_predictors = rcopy(R"(outputs$true_predictors)")
    intercept = rcopy(R"(params$intercept)")
    observations = rcopy(R"(params$observations)")


    # betas and basis
    beta_matrix = rcopy(R"(outputs$B)")
    basis_objs = rcopy(R"(outputs$basis_objs)")
    basis_values = rcopy(R"(outputs$basis_values)")
    time_domains    = rcopy(R"(params$time_domains)")

    # matrixes 
    X = rcopy(R"(outputs$X)")
    Y = rcopy(R"(outputs$Y)")
    Z = rcopy(R"(outputs$Z)")
    J = rcopy(R"(outputs$J)")
    W = rcopy(R"(outputs$W)");
    return X, Y, intercept, Z, J, W, beta_matrix, basis_objs, basis_values, predictors, true_predictors, time_domains
end