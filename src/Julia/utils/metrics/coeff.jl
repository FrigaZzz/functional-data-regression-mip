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