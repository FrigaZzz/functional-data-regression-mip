


# Modified run_model_and_save_outputs function
function run_model_and_save_outputs(model_name, simulation_name, setting_name, observations, observations_test, measurements, basis_functions, noise_snr, seed, λ, λ_group, M)
    # Include model file
    include(joinpath(project_root, "src", "Julia", "ols_vs_mip_models", model_name * ".jl"))

    # write the rput macro to pass the variable to R
    train_data = load_simulation_data(simulation_name, setting_name, project_root; observations = observations, measurements = measurements, basis_functions = basis_functions, noise_snr = noise_snr, seed =seed)

    Y_train, intercept, Z_train, beta_matrix_train, predictors_train, true_predictors_train = (
        train_data[:Y], train_data[:intercept], train_data[:Z], train_data[:B], train_data[:predictors], train_data[:true_predictors]
    )

    basis_values = train_data[:basis_values]

    test_data = load_simulation_data(simulation_name, setting_name, project_root; observations = observations_test, measurements = measurements, basis_functions = basis_functions, noise_snr = noise_snr, seed = seed + 1) 

    Y_test, Z_test = (
        test_data[:Y], test_data[:Z]
    )


    output_dir = define_output_dir(model_name, simulation_name, setting_name, observations, observations_test, basis_functions, measurements, λ, λ_group, M, noise_snr, seed, project_root)
    ensure_directory_exists(output_dir)


    beta_star = alpha = groups = nothing

    to_predict = sum(true_predictors_train)
    intercept = train_data[:intercept]

    BigM = ones(size(beta_matrix_train)) .*     100000   # or use beta_matrix_max_values
    BigM_ =  ones(size(beta_matrix_train)) .*  - 100000  # or use beta_matrix_min_values

    beta_star, alpha, groups = mip_functional_regression(Y_train, Z_train, BigM,BigM_; intercept = intercept != 0, group_limit= to_predict)

    # Compute Metrics and save to file
    performance_metrics = compute_metrics(Y_test, Z_test, beta_matrix_train, beta_star, alpha, groups, to_predict,basis_values)




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
    model_params."SNR" = [noise_snr]
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

    output_image_path = joinpath(output_dir, "curves")
    plots_path = joinpath(output_dir, "plots")

    plot_combined_predicted_curve(train_data[:beta_point_values], beta_star, train_data[:basis_values],train_data[:time_domains], output_image_path, false)
    generate_and_save_plots(train_data[:X], test_data[:X], Y_train, Y_test, train_data[:basis_values], train_data[:W], test_data[:W], Z_train, Z_test,plots_path)
    Y_pred = get_predictions(Z_test, beta_star, alpha)
    save_plots(Y_test, Y_pred, Z_test, beta_star, alpha, observations,plots_path)
end
project_root = dirname(@__DIR__)

include(joinpath(project_root, "setup", "init_env.jl"))

set_R_lib_path(project_root)
model_name = "l0.jl"
include(joinpath(project_root, "src", "Julia", "ols_vs_mip_models", model_name))
include(joinpath(project_root, "src",  "simulation.jl"))
include(joinpath(project_root, "src", "Julia", "utils", "file_management.jl"))
include(joinpath(project_root, "src", "Julia", "utils", "data_analysis.jl"))
include(joinpath(project_root, "src", "Julia", "utils", "plot.jl"))


if isempty(ARGS)
    # Predefined values
    run_model_and_save_outputs("l0", "3_predictors", "default", 2500, 500, 250, 6, [false,100], 1, 0.01, 0.001, 1000)

else
    # Values from command-line arguments
    run_model_and_save_outputs(parse_command_line_arguments(ARGS)...)
end