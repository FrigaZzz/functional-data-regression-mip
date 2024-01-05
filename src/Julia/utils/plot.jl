using Plots

"""
    plot_combined_predicted_curve(γ_matrix, γ_star, basis_values, time_domains, folder_path, show_plot=false)

Plot the combined predicted curve for γ_matrix and γ_star.

# Arguments
- `γ_matrix`: Coefficient matrix for the training data.
- `γ_star`: Coefficient matrix for the test data.
- `basis_values::Array`: Array of basis function values.
- `time_domains::Array`: Array of time domains.
- `folder_path::String`: Path to the folder where the plot will be saved.
- `show_plot::Bool=false`: Whether to display the plot.

# Returns
- `nothing`

"""
function plot_combined_predicted_curve(beta_point_values::Any, γ_star::Any, basis_values::Any, time_domains::Array, folder_path::Any, show_plot::Any=false;  upper_bound = nothing, lower_bound = nothing)
    mkpath(folder_path)  # Ensure the plots directory exists
    rm(folder_path, recursive=true)  # Clear the directory
    mkpath(folder_path)  # Recreate the directory
    predictors = size(γ_star, 1)

    # Calculate the number of rows for the layout based on the number of predictors
    n_rows = ceil(Int, predictors / 3)
    p = plot(layout = grid(n_rows, 3))

    for j in 1:predictors  # Loop over all curves
        curr_basis = basis_values[j,:,:]
        combined_curve_matrix = beta_point_values[j, :]
        combined_curve_star = curr_basis * γ_star[j, :]
        
        # Standardized time domain for plotting
        real_time = time_domains[j, :] 

        # Define subplot index
        subplot_index = j

        # Plot the combined predicted curve for γ_matrix and γ_star
        plot!(p[subplot_index], real_time, combined_curve_matrix, label="γ_matrix", color=:blue)
        plot!(p[subplot_index], real_time, combined_curve_star, label="γ_predicted", color=:red, linestyle=:dash)

        # Add the upper and lower bounds if they exist
        if !isnothing(upper_bound)
            plot!(p[subplot_index], real_time, curr_basis * upper_bound[j, :], label="Upper Bound", color=:green, linestyle=:dashdotdot)
        end
        if !isnothing(lower_bound)
            plot!(p[subplot_index], real_time, curr_basis * lower_bound[j, :], label="Lower Bound", color=:brown, linestyle=:dashdotdot)
        end
    end

    if show_plot
        display(p)
    end
    
    # Save the plot to a file
    savefig(p, joinpath(folder_path, "combined.png"))
end


# Function to plot multiple observations for a single predictor
function plot_predictor_observations(X, predictor, obs, plot_title)
    p = plot(X[1, predictor, :], label = "X Obs 1", title = plot_title)  # Start with the first observation
    for i in 2:obs
        plot!(p, X[i, predictor, :], label = "Obs $i")  # Add the rest
    end
    return p
end

# Function to plot the W matrix transformed by the basis values
function plot_transformed_W(basis_values, W, predictor, obs, plot_title)
    p = plot(basis_values[predictor,:,:] * W[1, predictor, :], label = "Obs 1", title = plot_title)  # Start with the first observation
    for i in 2:obs
        plot!(p, basis_values[predictor,:,:] * W[i, predictor, :], label = "Obs $i")  # Add the rest
    end
    return p
end

# Function to plot Z for multiple observations
function plot_Z_observations(Z, predictor, obs, plot_title)
    p = plot(Z[1, predictor, :], label = "Obs 1", title = plot_title)  # Start with the first observation
    for i in 2:obs
        plot!(p,  Z[i, predictor, :], label = "Obs $i")  # Add the rest
    end
    return p
end

# New function to generate and save plots
function generate_and_save_plots(X, X_test, Y, Y_test, basis_values, W, W_test, Z, Z_test, folder_path)
    mkpath(folder_path)  # Ensure the plots directory exists
    rm(folder_path, recursive=true)  # Clear the directory
    mkpath(folder_path)  # Recreate the directory
    

    predictor = 1
    obs = 20  # Adjust as needed

    # Create the plots
    pY = plot(Y, title = "Y Plot")
    pY_test = plot(Y_test, title = "Y Test Plot")
    p1 = plot_predictor_observations(X, predictor, obs, "Predictor Observations")
    p2 = plot_transformed_W(basis_values, W, predictor, obs, "Transformed W")
    p3 = plot_Z_observations(Z, predictor, obs, "Z Observations")
    p11 = plot_predictor_observations(X_test, predictor, obs, "Predictor Observations Test")
    p22 = plot_transformed_W(basis_values, W_test, predictor, obs, "Transformed W Test")
    p33 = plot_Z_observations(Z_test, predictor, obs, "Z Observations Test")

    # Combine plots side by side
    p_non_test = plot(pY, p1, p2, p3, layout = (1, 4), size = (1200, 300))
    p_test = plot(pY_test, p11, p22, p33, layout = (1, 4), size = (1200, 300))

    # Save plots to files
    file_path_non_test = joinpath(folder_path, "non_test.png")
    file_path_test = joinpath(folder_path, "test.png")
    savefig(p_non_test, file_path_non_test)
    savefig(p_test, file_path_test)
end



using Plots

function save_plots(Y_test, Y_pred, Z_test, beta_star, alpha_star, observations, output_dir)
    # 1. Scatter plot comparing true vs predicted values
    p1 = scatter(Y_test, Y_pred, xlabel="True Y", ylabel="Predicted Y", legend=false, title="True vs Predicted Y")
    plot!(Y_test, Y_test, color=:red)  # A y=x line for reference
    scatter_plot_path = joinpath(output_dir, "scatter_plot.png")
    savefig(p1, scatter_plot_path)

    # 2. Plot residuals
    residuals = Y_test - Y_pred
    p2 = scatter(1:observations, residuals, xlabel="Observation", ylabel="Residual", legend=false, title="Residuals")
    hline!([0], color=:red, label="Zero line")
    residuals_plot_path = joinpath(output_dir, "residuals_plot.png")
    savefig(p2, residuals_plot_path)
end

