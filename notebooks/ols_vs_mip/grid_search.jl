# Set up environment
# @__DIR__ is the directory of the current file
# We need to go up to the parent directory to find the project root
project_root = dirname(dirname(@__DIR__))
include(joinpath(project_root, "setup", "init_env.jl"))
set_R_lib_path(project_root)
using LinearAlgebra

# Load required packages
using Plots
using Statistics
using LinearAlgebra
using Random
include(joinpath(project_root, "src", "Julia", "utils", "data_analysis.jl"))
# Include simulation code
include(joinpath(project_root, "src", "simulation.jl"))

function run_lambda_test_grid(k_fold,observations, measurements, basis_functions, noise_snr, lambda,simulation_name,simulation_settings_file,lambda_values_to_test)

    include(joinpath(project_root, "src", "simulation.jl"))


    params_train = (
        observations=observations,
        measurements=measurements,
        basis_functions=basis_functions,
        noise_snr=noise_snr,
        seed=0
    )

    params_test = (
        observations=observations,
        measurements=measurements,
        basis_functions=basis_functions,
        noise_snr=noise_snr,
        seed=100
    )



    # Note: Use ... to unpack NamedTuple into keyword arguments
    output = load_simulation_data(simulation_name, simulation_settings_file, project_root; params_train...)

    output_test = load_simulation_data(simulation_name, simulation_settings_file, project_root; params_test...)


    # Grab the outputs from the R script

    predictors = Int(output[:predictors])
    true_predictors = output[:true_predictors]
    intercept = output[:intercept]
    observations = Int(output[:observations])

    # betas and basis
    beta_matrix = output[:B]
    basis_objs = output[:basis_objs]
    basis_values = output[:basis_values]
    time_domains = output[:time_domains]

    # matrixes 
    U = output[:U]
    X = output[:X]
    Y = output[:Y]
    Z = output[:Z]
    J = output[:J]
    W = output[:W]


    X_test = output_test[:X]
    Y_test = output_test[:Y]
    Z_test = output_test[:Z]
    J_test = output_test[:J]
    W_test = output_test[:W]
    beta_matrix_test = output_test[:B]



    model_name = "l0_and_l2"
    model_file_path = joinpath(project_root, "src", "Julia", "ols_vs_mip_models", model_name * ".jl")

    include(model_file_path)

    beta_matrix_max_values = maximum(beta_matrix, dims=2)
    beta_matrix_min_values = minimum(beta_matrix, dims=2)

    BigM = ones(size(beta_matrix)) .* 100000   # or use beta_matrix_max_values
    BigM_ = ones(size(beta_matrix)) .* -100000  # or use beta_matrix_min_values

    to_predict = sum(true_predictors)
    metrics = Dict()

    # Split data into k folds
    indices = collect(1:length(Y))
    Random.shuffle!(indices)
    folds = collect(Iterators.partition(indices, div(length(Y), k_fold)))

    avg_mse_per_lambda = Dict{Float64,Float64}()
    avg_ise_per_lambda = Dict{Float64,Float64}()

    for lambda in lambda_values_to_test
        println("lambda: ", lambda)
        fold_metrics = []
        fold_mse = []
        fold_ise = []
        for fold in 1:k_fold
            # Use fold for validation, rest for training
            validation_indices = folds[fold]
            train_indices = setdiff(indices, validation_indices)

            Y_train, Z_train = Y[train_indices], Z[train_indices, :, :]
            Y_val, Z_val = Y[validation_indices], Z[validation_indices, :, :]

            beta_star, alpha_star, groups = mip_functional_regression(Y_train, Z_train, BigM, BigM_;
                intercept=output[:intercept] != 0,
                group_limit=to_predict,
                lambda=lambda)

            performance_metrics = compute_metrics(Y_val, Z_val, beta_matrix_test, beta_star, alpha_star, groups, predictors)
            push!(fold_metrics, performance_metrics)
            push!(fold_ise, performance_metrics["ISE_Coefficients"])  # Assuming compute_metrics returns a dictionary with "MSE"

            push!(fold_mse, performance_metrics["MSE_Predictions"])  # Assuming compute_metrics returns a dictionary with "MSE"
        end
        # average the metrics over the folds that is an array of Dictioranies
        avg_metrics = Dict()
        for i in 1:length(fold_metrics)
            for (key, value) in fold_metrics[i]
                if haskey(avg_metrics, key)
                    avg_metrics[key] += value
                else
                    avg_metrics[key] = value
                end
            end
        end        
        # now put average values in the dictionary
        for (key, value) in avg_metrics
            avg_metrics[key] = value / length(fold_metrics)
        end
        metrics[lambda] = [avg_metrics]

        avg_mse = mean(fold_mse)
        avg_ise = mean(fold_ise)

        avg_mse_per_lambda[lambda] = avg_mse
        avg_ise_per_lambda[lambda] = avg_ise

    end
    print(avg_ise_per_lambda)
    # Choose lambda with lowest average MSE
    # Find the key corresponding to the minimum value in avg_mse_per_lambda
    best_lambda_mse_ix = argmin(collect(values(avg_mse_per_lambda)))
    best_lambda_mse = collect(keys(avg_mse_per_lambda))[best_lambda_mse_ix]
    # Find the key corresponding to the minimum value in avg_ise_per_lambda
    best_lambda_ise_ix = argmin(collect(values(avg_ise_per_lambda)))
    best_lambda_ise = collect(keys(avg_ise_per_lambda))[best_lambda_ise_ix]
    println("\n")

    println("Best selected lambda MSE: ", best_lambda_mse)
    println("Corresponding average MSE: ", avg_mse_per_lambda[best_lambda_mse])
    println("Corresponding average ISE: ", avg_ise_per_lambda[best_lambda_mse])

    println("Best lambda for ISE: ", best_lambda_ise)
    println("Corresponding average MSE: ", avg_mse_per_lambda[best_lambda_ise])
    println("Corresponding average ISE: ", avg_ise_per_lambda[best_lambda_ise])


    # Save results to a JSON file
    json_file_path = joinpath(project_root, "docs", "graphs", "sim3", "lambda", "results.json")
    open(json_file_path, "w") do io
        JSON.print(io, metrics)
    end
    return best_lambda_mse, best_lambda_ise, avg_mse_per_lambda[best_lambda_mse], avg_ise_per_lambda[best_lambda_mse], avg_mse_per_lambda[best_lambda_ise], avg_ise_per_lambda[best_lambda_ise]
end

# Example usage with 5-fold cross-validation
let
    simulation_name = "paper2"
    simulation_settings_file = "10_pred_sim2"
    k_folds = 10
    measurements = 50
    observations = 900
    basis_functions = 5
    noise_snr=  [100,100]
    lambda = 0.1
    lambda_values_to_test = range(0, 100, length=20)

    best_lambda_mse, best_lambda_ise, mse_best_lambda_mse, ise_best_lambda_mse, mse_best_lambda_ise, ise_best_lambda_ise = redirect_stdout(devnull) do
        run_lambda_test_grid(k_folds, observations, measurements, basis_functions, noise_snr, lambda, simulation_name, simulation_settings_file,lambda_values_to_test)
    end

    println("Best selected lambda MSE: ", best_lambda_mse)
    println("Corresponding average MSE: ", mse_best_lambda_mse)
    println("Corresponding average ISE: ", ise_best_lambda_mse)
    
    println("Best lambda for ISE: ", best_lambda_ise)
    println("Corresponding average MSE: ", mse_best_lambda_ise)
    println("Corresponding average ISE: ", ise_best_lambda_ise)
end;