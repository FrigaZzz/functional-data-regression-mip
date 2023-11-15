using Plots

function plot_combined_predicted_curve(γ_matrix, γ_star, basis_values, time_domain, folder_path, show_plot =false )
    mkpath(folder_path)  # Ensure the plots directory exists
    rm(folder_path, recursive=true)  # Clear the directory
    mkpath(folder_path)  # Recreate the directory
    
    for j in 1:size(γ_matrix, 1)  # Loop over all curves

        
        combined_curve_matrix = basis_values * γ_matrix[j, :]
        combined_curve_star = basis_values * γ_star[j, :]
        
        # Plot the combined predicted curve for γ_matrix and γ_star
        p = plot(time_domain, combined_curve_matrix, label="γ_matrix", color=:blue)
        plot!(time_domain, combined_curve_star, label="γ_predicted", color=:red, linestyle=:dash)
        
        if(show_plot)
            display(p)
        end
        # Save the plot to a file
        savefig(p, joinpath(folder_path, "predictor$(j)_combined_predicted_curve.png"))
    end
end


function plot_combined_predicted_curve_paper(γ_matrix, γ_star, basis_values,  time_domains, folder_path, show_plot=false)
    mkpath(folder_path)  # Ensure the plots directory exists
    rm(folder_path, recursive=true)  # Clear the directory
    mkpath(folder_path)  # Recreate the directory
    predictors =  size(γ_matrix, 1)
    for j in 1:predictors  # Loop over all curves
        curr_basis = basis_values[j,:,:]
        # reshape curr_basis
        combined_curve_matrix = curr_basis * γ_matrix[j, :]
        combined_curve_star = curr_basis * γ_star[j, :]
        
        
        # Plot the combined predicted curve for γ_matrix and γ_star
        p = plot(time_domains[[j]], combined_curve_matrix, label="γ_matrix", color=:blue)
        plot!(p,time_domains[[j]], combined_curve_star, label="γ_predicted", color=:red, linestyle=:dash)
        

        if(show_plot)
            display(p)
        end
        
        # Save the plot to a file
        savefig(p, joinpath(folder_path, "predictor$(j)_combined_predicted_curve.png"))
    end
end