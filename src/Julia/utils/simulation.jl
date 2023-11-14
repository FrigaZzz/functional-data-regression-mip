using RCall


# Function to load R simulation data
function load_simulation_data(simulation_name, project_root, observations, measurements, basis_functions, error_sd, seed)
    # clear R environment
    # load new simulation data
    simulation_file_path = joinpath(project_root, "simulations", "settings", simulation_name, "setting"  * ".R")
    @rput observations measurements basis_functions error_sd seed

    @rput simulation_file_path
    R"""source(simulation_file_path)"""
    # Load simulation data, setup directories, etc., now using the passed simulation_name
    # Setup variables from R environment
    predictors = rcopy(R"(params$predictors)")
    true_predictors = rcopy(R"(true_predictors)")
    intercept = rcopy(R"(params$intercept)")

    # betas and basis
    beta_matrix = rcopy(R"(B)")
    basis_obj = rcopy(R"(basis_obj)")
    basis_values = rcopy(R"(basis_values)")

    # matrixes 
    X = rcopy(R"(X)")
    Y = rcopy(R"(Y)")
    Z = rcopy(R"(Z)")
    J = rcopy(R"(J)")
    W = rcopy(R"(W)")
    return X, Y, intercept, Z, J, W, beta_matrix, basis_obj, basis_values, predictors, true_predictors
end