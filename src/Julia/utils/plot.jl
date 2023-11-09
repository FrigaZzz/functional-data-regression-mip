using Plots

function plot_combined_predicted_curve(γ_matrix, γ_star, basis_values, m, r, tb, folder_path )
    mkpath(folder_path)  # Ensure the plots directory exists
    rm(folder_path, recursive=true)  # Clear the directory
    mkpath(folder_path)  # Recreate the directory

    t_vals = LinRange(0, 1, tb)
    
    for j in 1:size(γ_matrix, 1)  # Loop over all curves
        p = plot(; legend=:topright, xlabel="t", ylabel="Combined Coefficient Value", title="Predictor $j: Combined Predicted Curve")
        
        combined_curve_matrix = basis_values * γ_matrix[j, :]
        combined_curve_star = basis_values * γ_star[j, :]
        
        # Plot the combined predicted curve for γ_matrix and γ_star
        plot!(t_vals, combined_curve_matrix, label="γ_matrix", color=:blue)
        plot!(t_vals, combined_curve_star, label="γ_predicted", color=:red, linestyle=:dash)
        
        # Save the plot to a file
        savefig(p, joinpath(folder_path, "predictor$(j)_combined_predicted_curve.png"))
    end
end