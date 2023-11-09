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

# Function to process results and generate plots
function process_and_plot(Y, Z, alpha_star, beta_matrix, beta_star, basis_values, n, p, r, tb)
    # Compute and display DataFrame
    df = compute_dataframe(beta_matrix, beta_star, p)

    # Save DataFrame to a CSV file
    df_file_path = joinpath(project_root, "outputs", "logs", simulation_name, "error_df.csv")
    CSV.write(df_file_path, df)

    # Compute evaluation metrics
    mse, rmse, mae, relative_error = compute_evaluation_metrics(beta_matrix, beta_star)

    # Save performance metrics to a text file
    performance_file_path = joinpath(project_root, "outputs", "logs", simulation_name, "performance.txt")
    open(performance_file_path, "w") do io
        write(io, "mse: $mse\nrmse: $rmse\nmae: $mae\nrelative_error: $relative_error\n")
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

# Function to compute evaluation metrics
function compute_evaluation_metrics(beta_matrix, beta_star)
    diff_matrix = beta_matrix - beta_star

    # Evaluation Metrics
    mse = sum((diff_matrix) .^ 2) / length(beta_star)  # Mean Squared Error
    rmse = sqrt(mse)                                  # Root Mean Squared Error
    mae = sum(abs.(diff_matrix)) / length(beta_star)  # Mean Absolute Error
    relative_error = norm(diff_matrix) / norm(beta_star)  # Relative Error

    return mse, rmse, mae, relative_error
end



# Function to generate and save plots
function generate_and_save_plots(Y, Y_pred, beta_matrix, beta_star, basis_values, n, p, r, tb)
    
    p1 = scatter(Y, Y_pred, xlabel="True Y", ylabel="Predicted Y", legend=false, title="True vs Predicted Y")
    # plot the Y as red line over the predicted Y
    plot!(p1,Y, Y, color=:red)  # A y=x line for reference
    save_plot(p1, "true_vs_predicted")

    residuals = Y - Y_pred
    p2 = scatter(1:n, residuals, xlabel="Observation", ylabel="Residual", legend=false, title="Residuals")
    hline!(p2, [0], color=:red, label="Zero line")
    save_plot(p2, "residuals")

    # Plot combined predicted curve
    plot_combined_predicted_curve(beta_matrix, beta_star, basis_values, p, r, tb, joinpath(project_root,"outputs", "plots", simulation_name))
end



# Run the full model process
run_model_and_save_outputs()
