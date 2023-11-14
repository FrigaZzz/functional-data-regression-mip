

project_root = dirname(@__DIR__)
using Pkg
Pkg.activate(project_root)
Pkg.instantiate()


using RCall

R"""
    .libPaths(c("C:/Users/User/AppData/Local/R/cache/R/renv/library/functional-data-regression-mip-e3349204/R-4.3/x86_64-w64-mingw32",
    .libPaths()))
"""


# Modified run_model_and_save_outputs function
function run_model_and_save_outputs(model_name, simulation_name, observations, observations_test, measurements,  basis_functions, error_sd, seed, λ, λ_group, M, override_true_predictors=nothing)
    # Include model file
    include(joinpath(project_root, "src", "Julia", "models", model_name))


    # Define your input variables
    # Setup R environment parameters using rput
    # write the rput macro to pass the variable to R
    X_train, Y_train, intercept, Z_train, J_train, W_train, beta_matrix_train, basis_obj_train, basis_values_train, predictors_train, true_predictors_train = load_simulation_data(simulation_name, project_root, observations, measurements, basis_functions, error_sd, seed)

    X_test, Y_test, _, Z_test, J_test, W_test, beta_matrix_test, basis_obj_test, basis_values_test, predictors_test, _ = load_simulation_data(simulation_name, project_root, observations_test, measurements, basis_functions, error_sd, seed)
    
    output_dir = define_output_dir(simulation_name, observations, observations_test, basis_functions, measurements, λ, λ_group, M, error_sd, seed, project_root)
    ensure_directory_exists(output_dir)


    beta_star = alpha = groups = nothing

    to_predict = sum(true_predictors_train) 
    # if(override_true_predictors != nothing)
    #     to_predict = override_true_predictors
    # end
    beta_star, alpha, groups = mip_functional_regression(Y_train, Z_train, λ, λ_group, M, to_predict)

    # Compute Metrics and save to file
    performance_metrics = compute_metrics(Y_test, Z_test, beta_matrix_train, beta_star, alpha, groups, true_predictors_train)
    save_performance_evaluation(output_dir, performance_metrics)

    # Save model parameters
    model_params = Dict(
        "Observations" => [observations],
        "Predictors" => [predictors_train],
        "BasisFunctions" => [basis_functions],
        "Measurements" => [measurements],
        "Lambda" => [λ],
        "LambdaGroup" => [λ_group],
        "M" => [M],
        "ErrorSD" => [error_sd],
        "Seed" => [seed]
    )
    save_model_parameters(output_dir, model_params)

    # Save model results
    model_results = Dict(
        "beta_star" => [beta_star],
        "True Intercept" => [[intercept]],
        "Evaluated Intercept" => [[alpha]],
        "groups" => [groups],
        "true_predictors" => [true_predictors_train],
        "beta_matrix" => [beta_matrix_train],
    )
    save_model_results(output_dir, model_results)


end

model_name = "l0_and_limit.jl"
include(joinpath(project_root, "src", "Julia", "models", model_name))
include(joinpath(project_root, "src", "Julia", "utils", "simulation.jl"))
include(joinpath(project_root, "src", "Julia", "utils", "file_management.jl"))
include(joinpath(project_root, "src", "Julia", "utils", "data_analysis.jl"))


if isempty(ARGS)
    # Predefined values
    override_true_predictors = 2 # else nothing
    run_model_and_save_outputs("l0_and_limit.jl", "paper", 500, 250, 100, 6, 0.001, 1, 0.01, 0.2, 1000,override_true_predictors)
    # run_model_and_save_outputs("l0_and_limit.jl", "5_predictors", 2500, 500, 250, 6, 0.001, 1, 0.01, 0.001, 1000)

else
    # Values from command-line arguments
    run_model_and_save_outputs(parse_command_line_arguments(ARGS))
end