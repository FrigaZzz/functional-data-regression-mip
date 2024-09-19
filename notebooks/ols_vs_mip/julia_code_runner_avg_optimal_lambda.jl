    # Set up environment
    # @__DIR__ is the directory of the current file
    # We need to go up to the parent directory to find the project root
    project_root = dirname(dirname(@__DIR__))
    # include(joinpath(project_root, "setup", "init_env.jl"))
    # set_R_lib_path(project_root)

    # Load required packages
    using Plots
    using Statistics
    # Include simulation code
    include(joinpath(project_root, "src", "simulation.jl"))
    # Load data analysis utilities
    using LinearAlgebra
    include(joinpath(project_root, "src", "Julia", "utils", "data_analysis.jl"))

    # Include plot file
    plot_file_path = joinpath(project_root, "src", "Julia","utils", "plot.jl")
    include(plot_file_path)   
    # Set output folder for plots
    
    # Include model file
    model_name = "l0_and_l2"
    model_file_path = joinpath(project_root, "src", "Julia",    "ols_vs_mip_models", model_name *".jl")
    include(model_file_path)
    # Define OLS solution function
    function ols_solution(Y, Z)
        Z_reshaped = reshape(Z, :, size(Z, 2) * size(Z, 3))
        Z_with_intercept = hcat(ones(size(Z_reshaped, 1)), Z_reshaped)
        beta_hat = (Z_with_intercept' * Z_with_intercept) \ (Z_with_intercept' * Y)
        beta_hat_matrix = reshape(beta_hat[2:end], size(Z, 2), size(Z, 3))
        return beta_hat_matrix
    end


    using Statistics

    function compute_auc(fpr, tpr)
        auc = 0
        for i in 1:length(fpr)-1
            auc += (fpr[i+1] - fpr[i]) * (tpr[i+1] + tpr[i]) / 2
        end
        return auc
    end

    # Function to accumulate metrics
    function accumulate_metrics(metrics_dict, cumulative_dict)
        for (key, value) in metrics_dict
            if typeof(value) == Float64
                cumulative_dict[key] = get(cumulative_dict, key, 0.0) + value
            end
        end
    end

    function compute_classification_scores_all_simulations(groups_list, true_predictors)
        recall_list, fpr_list, precision_list, f1_score_list, roc_auc_list = [], [], [], [], []
    
        for groups in groups_list
            # Compute True Positives (TP), False Positives (FP), True Negatives (TN), False Negatives (FN)
            TP = sum((true_predictors .== 1.0) .& (groups .== 1.0))
            FP = sum((true_predictors .== 0.0) .& (groups .== 1.0))
            TN = sum((true_predictors .== 0.0) .& (groups .== 0.0))
            FN = sum((true_predictors .== 1.0) .& (groups .== 0.0))
    
            # Compute Recall
            recall = TP / (TP + FN) # TPR
            fpr = FP / (FP + TN)

            # Compute Precision
            precision = TP / (TP + FP) #
    
            # Compute F1 Score
            f1_score = 2 * (precision * recall) / (precision + recall)
            auc = recall * (1 - fpr) # AUC for current simulation

            # Compute ROC AUC
            push!(fpr_list, fpr)
            push!(recall_list, recall)
            push!(precision_list, precision)
            push!(f1_score_list, f1_score)
            push!(roc_auc_list, auc)

        end

        return recall_list, precision_list, f1_score_list, roc_auc_list
    end
    



    # Define simulation parameters
    simulation_name = "paper2"
    simulation_settings_file = "default"
    output_folder = joinpath(project_root, "outputs", "plots",    simulation_name)  

    # Define the number of simulations
    N_simulations = 1
    basis_functions = 6
    measurements = 300
    observations = 500
    noise_snr=  [1,11]
    basis = range(1,30, length=30)
    base = 1
    coeff = 1
    lambda =base * 10.0^coeff



    # Initialize variables to store cumulative performance metrics
    # Initialize variables to store cumulative performance metrics
    cumulative_metrics_estimate = Dict{String, Float64}()
    cumulative_metrics_ols = Dict{String, Float64}()
    cumulative_metrics_real = Dict{String, Float64}()
    corrected_predictors = 0.0
    true_predictors = 
    groups_list =[]

    # Loop over the simulations
    for i in 1:N_simulations
        # i= i
        println("Simulation $i")
        global basis_functions, measurements, observations, noise_snr,true_predictors,groups_list,base,coeff

        println("LambdaLambdaLambdaLambdaLambdaLambda: ", lambda)
        params_train = (
            observations = observations,
            measurements = measurements,
            basis_functions = basis_functions,
            noise_snr = noise_snr,
            seed = i
            )

        params_test = (
            observations = observations,
            measurements = measurements,
            basis_functions = basis_functions,
            noise_snr = noise_snr,
            seed = i*10
        )

        # Load simulation data
        output = load_simulation_data(simulation_name, simulation_settings_file, project_root; params_train...)
        output_test = load_simulation_data(simulation_name, simulation_settings_file, project_root; params_test...)

        predictors = Int(output[:predictors])
        true_predictors = output[:true_predictors]
        intercept = output[:intercept]
        observations = Int(output[:observations])

        beta_matrix  = output[:B]
        basis_objs   = output[:basis_objs]
        basis_values = output[:basis_values]
        time_domains = output[:time_domains]

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
        beta_matrix_test  = output_test[:B]



        # Calculate maximum and minimum values of beta_matrix
        beta_matrix_max_values = maximum(beta_matrix, dims = 2)
        beta_matrix_min_values = minimum(beta_matrix, dims = 2)

        # Set BigM values
        BigM =ones(size(beta_matrix))  .* 30000   # or use   beta_matrix_max_values
        BigM_ =  ones(size(beta_matrix))  .* -30000  # or use    beta_matrix_min_values

        # Perform MIP functional regression
        to_predict = sum(true_predictors)
        if(model_name == "l0")
            beta_star, alpha_star, groups = mip_functional_regression(Y, Z, BigM,   BigM_; intercept = output[:intercept] != 0, group_limit = to_predict)
        else
            beta_star, alpha_star, groups = mip_functional_regression(Y, Z, BigM,   BigM_; intercept = output[:intercept] != 0, group_limit = to_predict, lambda = lambda)
        end

    
        # Calculate OLS coefficients on training data (omit code for brevity)
        beta_ols = ols_solution(Y, Z)
        print("Beta_ols: ", beta_ols)
        print("Beta_star: ", beta_star)
    
    # Compute performance metrics for testing data
    performance_metrics_real_test = compute_metrics(Y_test, Z_test, beta_matrix_test, beta_matrix, alpha_star, groups, predictors, basis_values)
    performance_metrics_estimate_test = compute_metrics(Y_test, Z_test, beta_matrix_test, beta_star, alpha_star, groups, predictors, basis_values)
    performance_metrics_ols_test = compute_metrics(Y_test, Z_test, beta_matrix_test, beta_ols, alpha_star, groups, predictors, basis_values)

    

    # Plot combined predicted curve
    beta_point_values = output_test[:beta_point_values]
    plot_combined_predicted_curve(beta_point_values, beta_star,   basis_values,   time_domains, output_folder, true)     
    # Update the number of corrected execution which means adding 1 if 
    # all the predictors are correctly estimated, which means that
    # groups = true_predictors position-wise, then +1


    global corrected_predictors
        if sum(groups .* true_predictors) == to_predict
            corrected_predictors += 1
        end
    # add groups into groups_list 
    push!(groups_list, groups)

    #Accumulate performance metrics
    accumulate_metrics(performance_metrics_real_test, cumulative_metrics_real)
    accumulate_metrics(performance_metrics_estimate_test, cumulative_metrics_estimate)
    accumulate_metrics(performance_metrics_ols_test, cumulative_metrics_ols)
end

    # Average the cumulative performance metrics over all simulations
    average_metrics_estimate = Dict{String, Float64}()
    average_metrics_ols = Dict{String, Float64}()
    average_metrics_real = Dict{String, Float64}()
    for (key, value) in cumulative_metrics_estimate
        average_metrics_estimate[key] = value / N_simulations
    end

    for (key, value) in cumulative_metrics_ols
        average_metrics_ols[key] = value / N_simulations
    end

    for (key, value) in cumulative_metrics_real
        average_metrics_real[key] = value / N_simulations
    end
    # Print or analyze the average performance metrics as needed

    # Print or analyze the average MSE as needed
    println("Average MSE for the real model: ", average_metrics_real)

    println("Average MSE for the estimated model: ", average_metrics_estimate)
    println("Average MSE for the OLS model: ", average_metrics_ols)
    println("Corrected predictors in : ", corrected_predictors, " execitions out of ", N_simulations )
    recall_list, precision_list, f1_score_list, roc_auc_list = compute_classification_scores_all_simulations(groups_list, true_predictors)

    # Aggregate results over all simulations
    mean_recall = mean(recall_list)
    mean_precision = mean(precision_list)
    mean_f1_score = mean(f1_score_list)
    mean_roc_auc = mean(roc_auc_list)
    
    println("Mean Recall: ", mean_recall)
    println("Mean Precision: ", mean_precision)
    println("Mean F1 Score: ", mean_f1_score)
    # println("Mean ROC AUC: ", mean_roc_auc)
    #clear all the created variables in Julia
