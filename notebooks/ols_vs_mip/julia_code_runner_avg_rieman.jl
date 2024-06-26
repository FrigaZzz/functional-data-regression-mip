# Set up environment
# @__DIR__ is the directory of the current file
# We need to go up to the parent directory to find the project root
project_root = dirname(dirname(@__DIR__))
include(joinpath(project_root, "setup", "init_env.jl"))
set_R_lib_path(project_root)

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
model_name = "l0"
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

# Define simulation parameters
simulation_name = "paper"
simulation_settings_file = "default"
output_folder = joinpath(project_root, "outputs", "plots",    simulation_name)  
measurements = 50
basis_functions = 6

params_train = (
    observations = 50,
    measurements = measurements,
    basis_functions = basis_functions,
    noise_snr = [100,1000],
    seed = 100
)

params_test = (
    observations = 50,
    measurements = measurements,
    basis_functions = basis_functions,
    noise_snr = [100,1000],
    seed = 300
)

# Define the number of simulations
N_simulations = 1

# Initialize variables to store cumulative performance metrics
# Initialize variables to store cumulative performance metrics
cumulative_metrics_estimate = Dict{String, Float64}()
cumulative_metrics_ols = Dict{String, Float64}()
cumulative_metrics_real = Dict{String, Float64}()
corrected_predictors = 0.0

# Loop over the simulations
for i in 1:N_simulations
    println("Simulation $i")

    params_train = (
    observations = 300,
    measurements = measurements,
    basis_functions = basis_functions,
    noise_snr = [100,1000],
    seed = i
    )

    params_test = (
        observations = 100,
        measurements = measurements,
        basis_functions = basis_functions,
        noise_snr = [100,1000],
        seed = i*10
    )

    # Load simulation data
    output = load_simulation_reiman(simulation_name, simulation_settings_file, project_root; params_train...)
    output_test = load_simulation_reiman(simulation_name, simulation_settings_file, project_root; params_test...)

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
    BigM =ones(size(beta_matrix))  .* 30000000   # or use   beta_matrix_max_values
    BigM_ =  ones(size(beta_matrix))  .* -30000000  # or use    beta_matrix_min_values

    # Perform MIP functional regression
    to_predict = sum(true_predictors)
    beta_star, alpha_star, groups = mip_functional_regression(Y, Z, BigM,   BigM_; intercept = output[:intercept] != 0, group_limit = to_predict)

    # Calculate OLS coefficients on training data (omit code for brevity)
    beta_ols = ols_solution(Y, Z)

   
   # Compute performance metrics for testing data
   performance_metrics_real_test = compute_metrics(Y_test, Z_test, beta_matrix_test, beta_matrix, alpha_star, groups, predictors)
   performance_metrics_estimate_test = compute_metrics(Y_test, Z_test, beta_matrix_test, beta_star, alpha_star, groups, predictors)
   performance_metrics_ols_test = compute_metrics(Y_test, Z_test, beta_matrix_test, beta_ols, alpha_star, groups, predictors)

   

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

   #Accumulate performance metrics
   for (key, value) in performance_metrics_real_test
    cumulative_metrics_real[key] = get(cumulative_metrics_real, key, 0.0) + value
end

   for (key, value) in performance_metrics_estimate_test
       cumulative_metrics_estimate[key] = get(cumulative_metrics_estimate, key, 0.0) + value
   end

   for (key, value) in performance_metrics_ols_test
       cumulative_metrics_ols[key] = get(cumulative_metrics_ols, key, 0.0) + value
   end
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
print("Corrected predictors in : ", corrected_predictors, " execitions out of ", N_simulations )