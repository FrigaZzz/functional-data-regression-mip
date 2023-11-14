using Statistics


# Metrics for Coefficient Evaluation
function mse_coefficients(beta_matrix, beta_star)
    mean((beta_matrix - beta_star) .^ 2)
end

function rmse_coefficients(beta_matrix, beta_star)
    sqrt(mse_coefficients(beta_matrix, beta_star))
end