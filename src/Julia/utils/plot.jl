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
function plot_combined_predicted_curve(beta_point_values::Any, γ_star::Any, basis_values::Any, time_domains::Array, folder_path::Any, show_plot::Any=false)

    mkpath(folder_path)  # Ensure the plots directory exists
    rm(folder_path, recursive=true)  # Clear the directory
    mkpath(folder_path)  # Recreate the directory
    predictors =  size(γ_star, 1)
    for j in 1:predictors  # Loop over all curves
        curr_basis = basis_values[j,:,:]
        # reshape curr_basis
        combined_curve_matrix = beta_point_values[j, :]
        combined_curve_star = curr_basis * γ_star[j, :]
        
        
        # Plot the combined predicted curve for γ_matrix and γ_star
        p = plot(time_domains[[j]], combined_curve_matrix, label="γ_matrix", color=:blue)
        plot!(p,time_domains[[j]], combined_curve_star, label="γ_predicted", color=:red, linestyle=:dash)
        

        if(show_plot)
            display(p)
        end
        
        # Save the plot to a file
        savefig(p, joinpath(folder_path, "predictor$(j).png"))
    end
end