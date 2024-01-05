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
include(joinpath(project_root, "simulations", "simulation.jl"))

# Define simulation parameters
simulation_name = "3_predictors"
simulation_settings_file = "default"

measurements = 100
basis_functions = 6

params_train = (
    observations = 150,
    measurements = measurements,
    basis_functions = basis_functions,
    noise_snr = [100,1000],
    seed = 1
)

params_test = (
    observations = 100,
    measurements = measurements,
    basis_functions = basis_functions,
    noise_snr = [100,1000],
    seed = 300
)

# Load simulation data
output = load_simulation_data(simulation_name, simulation_settings_file, project_root; params_train...)
output_test = load_simulation_data(simulation_name, simulation_settings_file, project_root; params_test...)

# Extract outputs from the simulation
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

# Plot observed X and basis values
plot(X[1,1,:], label="OBSERVED X", legend=:topleft)
plot!(basis_values[1,:,:] * W[1, 1, :], label="W", legend=:topleft)

# Include model file
model_name = "simple_regressor"
model_file_path = joinpath(project_root, "src", "Julia","ols_vs_mip_models", model_name *".jl")
include(model_file_path)

# Calculate maximum and minimum values of beta_matrix
beta_matrix_max_values = maximum(beta_matrix, dims = 2)
beta_matrix_min_values = minimum(beta_matrix, dims = 2)

# Set BigM values
BigM = ones(predictors) .* 30000000   # or use beta_matrix_max_values
BigM_ =  ones(predictors) .* -30000000  # or use beta_matrix_min_values

# Perform MIP functional regression
to_predict = sum(true_predictors)
beta_star, alpha_star, groups = mip_functional_regression(Y, Z, BigM, BigM_; intercept = output[:intercept] != 0, group_limit = to_predict)

# Define OLS solution function
function ols_solution(Y, Z)
    Z_reshaped = reshape(Z, :, size(Z, 2) * size(Z, 3))
    Z_with_intercept = hcat(ones(size(Z_reshaped, 1)), Z_reshaped)
    beta_hat = (Z_with_intercept' * Z_with_intercept) \ (Z_with_intercept' * Y)
    beta_hat_matrix = reshape(beta_hat[2:end], size(Z, 2), size(Z, 3))
    return beta_hat_matrix
end

# Calculate OLS coefficients
beta_ols = ols_solution(Y, Z)

# Include plot file
plot_file_path = joinpath(project_root, "src", "Julia","utils", "plot.jl")
include(plot_file_path)

# Set output folder for plots
output_folder = joinpath(project_root, "outputs", "plots", simulation_name)

# Plot combined predicted curve
beta_point_values = output[:beta_point_values]
plot_combined_predicted_curve(beta_point_values, beta_ols, basis_values, time_domains, output_folder, true)

# Load data analysis utilities
using LinearAlgebra
include(joinpath(project_root, "src", "Julia", "utils", "data_analysis.jl"))

# Compute performance metrics
performance_metrics_real = compute_metrics(Y, Z, beta_matrix, beta_matrix, alpha_star, groups, predictors)
performance_metrics_estimate = compute_metrics(Y, Z, beta_matrix, beta_star, alpha_star, groups, predictors)
performance_metrics_ols = compute_metrics(Y, Z, beta_matrix, beta_ols, alpha_star, groups, predictors)
