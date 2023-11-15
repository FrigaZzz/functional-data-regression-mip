
# Modified run_model_and_save_outputs function
function run_model_and_save_outputs(model_name, simulation_name, setting_name, observations, observations_test, measurements, basis_functions, error_sd, seed, λ, λ_group, M, override_true_predictors=nothing)
    # Include model file
    include(joinpath(project_root, "src", "Julia", "models", model_name * ".jl"))

    # write the rput macro to pass the variable to R
    train_data = load_simulation_data(simulation_name, setting_name, project_root, observations, measurements, basis_functions, error_sd, seed)

    Y_train, intercept, Z_train, beta_matrix_train, predictors_train, true_predictors_train = (
        train_data[:Y], train_data[:intercept], train_data[:Z], train_data[:B], train_data[:predictors], train_data[:true_predictors]
    )

    test_data = load_simulation_data(simulation_name, setting_name, project_root, observations_test, measurements, basis_functions, error_sd, seed)

    Y_test, Z_test = (
        test_data[:Y], test_data[:Z]
    )


    output_dir = define_output_dir(model_name, simulation_name, observations, observations_test, basis_functions, measurements, λ, λ_group, M, error_sd, seed, project_root)
    ensure_directory_exists(output_dir)


    beta_star = alpha = groups = nothing

    to_predict = sum(true_predictors_train)

    beta_star, alpha, groups = mip_functional_regression(Y_train, Z_train, λ, λ_group, M, to_predict)

    # Compute Metrics and save to file
    performance_metrics = compute_metrics(Y_test, Z_test, beta_matrix_train, beta_star, alpha, groups, true_predictors_train)
    save_performance_evaluation(output_dir, performance_metrics)

    # Save model parameters
    # Initialize an empty DataFrame
    model_params = DataFrame()
    # extract only the last part of the path output_dir
    last = split(output_dir, "/")[end]
    # Add columns one by one in the desired order
    model_params."SimulationName" = [simulation_name]
    model_params."Directory" = [last]
    model_params."SettingName" = [setting_name]
    model_params."ModelName" = [model_name]
    model_params."Observations" = [observations]
    model_params."Predictors" = [predictors_train]
    model_params."BasisFunctions" = [basis_functions]
    model_params."Measurements" = [measurements]
    model_params."Lambda" = [λ]
    model_params."LambdaGroup" = [λ_group]
    model_params."M" = [M]
    model_params."ErrorSD" = [error_sd]
    model_params."Seed" = [seed]


    save_model_parameters(output_dir, model_params)

    # Save model results
    model_results = Dict(
        "beta_star" => [beta_star],
        "True Intercept" => [[intercept]],
        "Evaluated Intercept" => [[alpha]],
        "groups" => [groups],
        "true_predictors" => [true_predictors_train],
        "beta_matrix" => [beta_matrix_train]
    )
    save_model_results(output_dir, model_results)


end

include(joinpath(project_root, "setup", "init_env.jl"))

project_root = dirname(@__DIR__)
set_R_lib_path(project_root)
# model_name = "l0_and_limit.jl"
# include(joinpath(project_root, "src", "Julia", "models", model_name))
include(joinpath(project_root, "src", "Julia", "utils", "simulation.jl"))
include(joinpath(project_root, "src", "Julia", "utils", "file_management.jl"))
include(joinpath(project_root, "src", "Julia", "utils", "data_analysis.jl"))


if isempty(ARGS)
    # Predefined values
    run_model_and_save_outputs("l0_and_limit", "5_predictors", "default", 2500, 500, 250, 6, 0.01, 1, 0.01, 0.001, 1000, nothing)

else
    # Values from command-line arguments
    run_model_and_save_outputs(parse_command_line_arguments(ARGS))
end