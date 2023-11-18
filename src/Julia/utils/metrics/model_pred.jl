"""
    rmse_predictions(Y_true, Y_pred)

Compute the root mean squared error (RMSE) between the true values `Y_true` and the predicted values `Y_pred`.

# Arguments
- `Y_true::Vector`: A vector of true values.
- `Y_pred::Vector`: A vector of predicted values.

# Returns
- `rmse::Float64`: The RMSE between `Y_true` and `Y_pred`.
"""
function rmse_predictions(Y_true, Y_pred)
    n = length(Y_true)
    return sqrt(sum((Y_true[i] - Y_pred[i])^2 for i in 1:n) / n)
end

"""
    mae_predictions(Y_true, Y_pred)

Compute the mean absolute error (MAE) between the true values `Y_true` and the predicted values `Y_pred`.

# Arguments
- `Y_true::Vector`: A vector of true values.
- `Y_pred::Vector`: A vector of predicted values.

# Returns
- `mae::Float64`: The MAE between `Y_true` and `Y_pred`.
"""
function mae_predictions(Y_true, Y_pred)
    n = length(Y_true)
    mae = sum(abs(Y_true[i] - Y_pred[i]) for i in 1:n) / n
    return  mae
end
"""
    mse_predictions(Y_true, Y_pred)

Compute the mean squared error (MSE) between the true values `Y_true` and the predicted values `Y_pred`.

# Arguments
- `Y_true::Vector`: a vector of true values.
- `Y_pred::Vector`: a vector of predicted values.

# Returns
- `mse::Float64`: the mean squared error between `Y_true` and `Y_pred`.
"""
function mse_predictions(Y_true, Y_pred)
    n = length(Y_true)
    mse = sum((Y_true[i] - Y_pred[i])^2 for i in 1:n) / n
    return  mse
end

"""
    r_squared(Y_test, Y_pred)

Compute the R-squared value between the true values `Y_test` and the predicted values `Y_pred`.

# Arguments
- `Y_test::Vector`: a vector of true values.
- `Y_pred::Vector`: a vector of predicted values.

# Returns
- `r2::Float64`: the R-squared value between `Y_test` and `Y_pred`.
"""
function r_squared(Y_test, Y_pred)
    ss_res = sum((Y_test .- Y_pred) .^ 2)
    ss_tot = sum((Y_test .- mean(Y_test)) .^ 2)
    1 - ss_res / ss_tot
end

"""
    adjusted_r_squared(Y_true, Y_pred, p)

Compute the adjusted R-squared value between the true values `Y_true` and the predicted values `Y_pred`.

# Arguments
- `Y_true::Vector`: a vector of true values.
- `Y_pred::Vector`: a vector of predicted values.
- `p::Int`: the number of predictors in the model.

# Returns
- `adj_r2::Float64`: the adjusted R-squared value between `Y_true` and `Y_pred`.
"""
function adjusted_r_squared(Y_true, Y_pred, number_of_predictors)
    
    # Assuming y_actual and y_predicted are arrays of actual and    predicted values
    n = length(Y_true) # Number of observations
    p = number_of_predictors # Replace with the actual number of    predictors in your model

    # Calculate R-squared
    SS_res = sum((Y_true - Y_pred) .^ 2)
    SS_tot = sum((Y_true .- mean(Y_pred)) .^ 2)
    r_squared = 1 - (SS_res / SS_tot)

    # Calculate Adjusted R-squared
    adjusted_r_squared = 1 - ((1 - r_squared) * ((n - 1) / (n - p - 1)))
    return adjusted_r_squared
end

