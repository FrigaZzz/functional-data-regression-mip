using Statistics
"""
    mse_coefficients(beta_matrix, beta_star)

Compute the mean squared error between two coefficient matrices `beta_matrix` and `beta_star`.

# Arguments
- `beta_matrix::Matrix`: the estimated coefficient matrix.
- `beta_star::Matrix`: the true coefficient matrix.

# Returns
- `mse::Float64`: the mean squared error between `beta_matrix` and `beta_star`.
"""


# Metrics for Coefficient Evaluation
function mse_coefficients(beta_matrix, beta_star)
    mean((beta_matrix - beta_star) .^ 2)
end

"""
    rmse_coefficients(beta_matrix, beta_star)

Compute the root mean squared error (RMSE) between two coefficient matrices `beta_matrix` and `beta_star`.

# Arguments
- `beta_matrix::Matrix`: the estimated coefficient matrix.
- `beta_star::Matrix`: the true coefficient matrix.

# Returns
- `rmse::Float64`: the RMSE between `beta_matrix` and `beta_star`.
"""
function rmse_coefficients(beta_matrix, beta_star)
    sqrt(mse_coefficients(beta_matrix, beta_star))
end


function se_coefficients(beta_matrix, beta_star, basis_values) 
    predictors = size(basis_values, 1)
    se = zeros(predictors)
    for j in 1:predictors
        curr_basis = basis_values[j,:,:]
        combined_curve_matrix = curr_basis * beta_matrix[j, :]
        combined_curve_star = curr_basis * beta_star[j, :]

        # Compute the root mean square error for the jth predictor
        se[j] = sqrt(sum((combined_curve_matrix .- combined_curve_star) .^ 2))
    end
    return se
end