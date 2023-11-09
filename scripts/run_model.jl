using RCall
using DataFrames
using Statistics
using LinearAlgebra
using Plots
using CSV

# Constants
project_root = dirname(@__DIR__)
simulation_name = "functional_data_analysis_3_pred_monovariate"

# Include model file
model_name = "l0_and_limit_refined.jl"
include(joinpath(project_root, "src", "Julia", "models", model_name))

# Include plot utility file
plot_file_path = joinpath(project_root, "src", "Julia", "utils", "plot.jl")
include(plot_file_path)

# Function to load R simulation data
function load_simulation_data()
    simulation_file_path = joinpath(project_root, "simulations", simulation_name * ".R")
    @rput simulation_file_path
    R"""source(simulation_file_path)"""
end

# Function to create directory if it doesn't exist
function ensure_directory_exists(path::String)
    if !isdir(path)
        mkpath(path)
    end
end

# Function to save plots to a specified subdirectory
function save_plot(plt, plot_type::String)
    plot_dir = joinpath(project_root, "outputs", "plots", simulation_name)
    plot_path = joinpath(plot_dir, "$(plot_type).png")
    savefig(plt, plot_path)
end

# Function to run the model and save outputs
function run_model_and_save_outputs()
    load_simulation_data() # Load R data

    # Setup output directory for logs
    output_dir = joinpath(project_root, "outputs", "logs", simulation_name)
    ensure_directory_exists(output_dir)

    # Setup variables from R environment
    n = Int(rcopy(R"params$observations"))
    p = Int(rcopy(R"params$predictors"))
    r = Int(rcopy(R"params$basis_functions"))
    t = Int(rcopy(R"params$measurements"))
    tb =Int(rcopy(R"params$basis_measurements"))
    
    # betas and basis
    beta_matrix = rcopy(R"B")
    basis_obj = rcopy(R"basis_obj")
    basis_values = rcopy(R"basis_values")
    
    # matrixes 
    Y = rcopy(R"Y")
    Z = rcopy(R"Z")
    J = rcopy(R"J")
    W = rcopy(R"W")

    # Model parameters
    lambda = 0.1
    lambda_group = 10
    M = 1000

    
    # Create the directory if it doesn't exist
    if !isdir(output_dir)
        mkpath(output_dir)
    end

    # Define the file paths where the output will be saved
    output_file_path = joinpath(output_dir,  "results.txt")
    outputs_file_path = joinpath(output_dir, "outputs.txt")

    beta_star = alpha_star = group = nothing

    # Open the file for writing
    open(output_file_path, "w") do io
        # Redirect stdout to the file
        redirect_stdout(io) do
            # Run the model code
            beta_star, alpha_star, group = mip_functional_regression(Y, Z, lambda,lambda_group, M)
        end
    end

    # Open the file for writing
    open(outputs_file_path, "w") do io
        # Write the outputs to the file
        write(io, "beta_star: $beta_star\nalpha_star: $alpha_star\ngroup: $group\n")
    end

    # Continue with other processing and plots
    process_and_plot(Y, Z, alpha_star, beta_matrix, beta_star, basis_values, n, p, r, tb)
end

# Function to obtain predictions using the estimated coefficients
function get_predictions(X_coeff, gamma_star, alpha_star)
    n, p, r = size(X_coeff)
    _, r_star = size(gamma_star) # get the actual size of gamma_star

    predictions = zeros(n)
    for i in 1:n
        prediction = 0
        for j in 1:p
            for k in 1:min(r, r_star)  # use the smaller of r and r_star
                prediction += X_coeff[i, j, k] * gamma_star[j, k]
            end
        end
        predictions[i] = prediction
    end

    return predictions
end

using Statistics

# RMSE
function calculate_rmse(Y_true, Y_pred)
    n = length(Y_true)
    return sqrt(sum((Y_true[i] - Y_pred[i])^2 for i in 1:n) / n)
end

# MAE
function calculate_mae(Y_true, Y_pred)
    n = length(Y_true)
    mae = sum(abs(Y_true[i] - Y_pred[i]) for i in 1:n) / n
    return  mae
end

# MSE
function calculate_mse(Y_true, Y_pred)
    n = length(Y_true)
    mse = sum((Y_true[i] - Y_pred[i])^2 for i in 1:n) / n
    return  mse
end


# MAPE
function calculate_mape(Y_true, Y_pred)
    n = length(Y_true)
    return 100 * sum(abs((Y_true[i] - Y_pred[i]) / Y_true[i]) for i in 1:n if Y_true[i] != 0) / n
end

# MAD
function calculate_mad(Y_true, Y_pred)
    return median(abs.(Y_true - Y_pred))
end

# Adjusted R^2
function calculate_adjusted_r2(Y_true, Y_pred, p)
    n = length(Y_true)
    mse = sum((Y_true[i] - Y_pred[i])^2 for i in 1:n) / n
    var_y_true = var(Y_true)
    r2 = 1 - mse / var_y_true
    return 1 - (1 - r2) * (n - 1) / (n - p - 1)
end

# NRMSE
function calculate_nrmse(Y_true, Y_pred)
    rmse = calculate_rmse(Y_true, Y_pred)
    return rmse / (maximum(Y_true) - minimum(Y_true))
end

# ABC
function calculate_abc(Y_true, Y_pred)
    return sum(abs.(Y_true - Y_pred))
end


# Add a qualitative valuation
function valuate_metric(metric_value, good_threshold, excellent_threshold=nothing)
    if isnothing(excellent_threshold)
        return metric_value < good_threshold ? "Good" : "Poor"
    else
        if metric_value < excellent_threshold
            return "Excellent"
        elseif metric_value < good_threshold
            return "Good"
        else
            return "Poor"
        end
    end
end


# Function to process results and generate plots
function process_and_plot(Y, Z, alpha_star, beta_matrix, beta_star, basis_values, n, p, r, tb)
    # Compute and display DataFrame
    df = compute_dataframe(beta_matrix, beta_star, p)

    # Save DataFrame to a CSV file
    df_file_path = joinpath(project_root, "outputs", "logs", simulation_name, "error_df.csv")
    CSV.write(df_file_path, df)

    
    # Get predictions
    Y_pred = get_predictions(Z, beta_star, alpha_star)
    
    
    
    # Evaluate the model using new metrics
    rmse = calculate_rmse(Y, Y_pred)
    mape = calculate_mape(Y, Y_pred)
    mae = calculate_mae(Y, Y_pred)
    mse = calculate_mse(Y, Y_pred)
    mad = calculate_mad(Y, Y_pred)
    adjusted_r2 = calculate_adjusted_r2(Y, Y_pred, p)
    nrmse = calculate_nrmse(Y, Y_pred)
    abc = calculate_abc(Y, Y_pred)
    
    # Valuation of metrics
    println("MSE: $mse, MAE: $mae")
    println("RMSE: $rmse (Lower is better. Below 0.1 is excellent, but depends on data scale.)")
    println("MAPE: $mape% (Lower is better. Typically < 10% is good.)")
    println("MAD: $mad (Lower is better. Less sensitive to outliers than RMSE.)")
    println("Adjusted R2: $adjusted_r2 (Higher is better. Closer to 1 indicates a better fit.)")
    println("NRMSE: $nrmse (Lower is better. Below 0.1 is excellent, but depends on data scale.)")
    println("ABC: $abc (Lower is better. No standard threshold; context-dependent.)")
    

    # Example of using the valuation function
    rmse_valuation = valuate_metric(rmse, 0.2, 0.1)
    println("RMSE Valuation: $rmse_valuation")
    

    # Save performance metrics to a text file
    performance_file_path = joinpath(project_root, "outputs", "logs", simulation_name, "performance.txt")
    open(performance_file_path, "w") do io
        write(io, "mse: $mse\nrmse: $rmse\nmae: $mae\nmape: $mape%\nmad: $mad\nadjusted_r2: $adjusted_r2\nnrmse: $nrmse\nabc: $abc\n")
    end

    Y_pred = get_predictions( Z, beta_star, alpha_star)

    # Generate and save plots
    generate_and_save_plots(Y, Y_pred, beta_matrix, beta_star, basis_values, n, p, r, tb)
end

# Function to compute DataFrame for results
function compute_dataframe(beta_matrix, beta_star, p)
    df = DataFrame()

    for i in 1:size(beta_matrix, 2)
        beta_col = Symbol("beta_$i")
        beta_star_col = Symbol("beta_star_$i")
        percent_error_col = Symbol("%_err_$i")
        df[!, beta_col] = beta_matrix[:, i]
        df[!, beta_star_col] = beta_star[:, i]

        # Calculate the percentage error while handling division by zero, using NaN to represent undefined values
        df[!, percent_error_col] = ifelse.(
            df[!, beta_col] .== 0,
            0, # Could use `missing` or another value if NaN is not suitable
            round.(((df[!, beta_star_col] .- df[!, beta_col]) ./ df[!, beta_col]) .* 100, digits=2)
        )
    end

    return df
end



# Function to generate and save plots
function generate_and_save_plots(Y, Y_pred, beta_matrix, beta_star, basis_values, n, p, r, tb)
    
 

    # Plot combined predicted curve
    plot_combined_predicted_curve(beta_matrix, beta_star, basis_values, p, r, tb, joinpath(project_root,"outputs", "plots", simulation_name))

    p1 = scatter(Y, Y_pred, xlabel="True Y", ylabel="Predicted Y", legend=false, title="True vs Predicted Y")
    # plot the Y as red line over the predicted Y
    plot!(p1,Y, Y, color=:red)  # A y=x line for reference
    save_plot(p1, "true_vs_predicted")

    residuals = Y - Y_pred
    p2 = scatter(1:n, residuals, xlabel="Observation", ylabel="Residual", legend=false, title="Residuals")
    hline!(p2, [0], color=:red, label="Zero line")
    save_plot(p2, "residuals")
end



# Run the full model process
run_model_and_save_outputs()
