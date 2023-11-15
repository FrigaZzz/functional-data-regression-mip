using DataFrames
using CSV
using RCall
using LinearAlgebra
using Statistics

include(joinpath(project_root, "src", "Julia", "utils", "metrics", "metrics.jl"))

# Function to obtain predictions using the estimated coefficients
"""
    get_predictions(Z_coeff, beta_star, alpha_star)

Compute the predictions for a functional data regression model.

# Arguments
- `Z_coeff::Array`: an array of shape `(n, p, r)` containing the functional covariates.
- `beta_star::Array`: an array of shape `(p, r_star)` containing the estimated regression coefficients.
- `alpha_star::Float64`: the estimated intercept.

# Returns
- `predictions::Array`: an array of shape `(n,)` containing the predicted values.
"""
function get_predictions(Z_coeff, beta_star, alpha_star)
    n, p, r = size(Z_coeff)
    _, r_star = size(beta_star) # get the actual size of beta_star

    predictions = zeros(n)
    for i in 1:n
        prediction = 0
        for j in 1:p
            for k in 1:min(r, r_star)  # use the smaller of r and r_star
                prediction += Z_coeff[i, j, k] * beta_star[j, k]
            end
        end
        predictions[i] = prediction + alpha_star
    end

    return predictions
end

"""
    compute_metrics(Y_test, Z_coeff, beta_matrix, beta_star, alpha_star, groups, true_predictors_train)

Compute various evaluation metrics for a functional data regression model.

# Arguments
- `Y_test::Array`: Test response variable.
- `Z_coeff::Array`: Coefficients of the functional predictors.
- `beta_matrix::Array`: Matrix of coefficients for the functional predictors.
- `beta_star::Array`: True coefficients for the functional predictors.
- `alpha_star::Float64`: True intercept.
- `groups::Array`: Group labels for the functional predictors.
- `true_predictors_train::Array`: Matrix of true predictors.

# Returns
- `Dict`: A dictionary of evaluation metrics.

# Evaluation Metrics
- `MSE_Coefficients`: Mean Squared Error of the coefficients.
- `RMSE_Coefficients`: Root Mean Squared Error of the coefficients.
- `ISE_Coefficients`: Integrated Squared Error of the coefficients.
- `MSE_Predictions`: Mean Squared Error of the predictions.
- `RMSE_Predictions`: Root Mean Squared Error of the predictions.
- `MAE_Predictions`: Mean Absolute Error of the predictions.
- `R_squared`: R-squared value of the predictions.
- `Adjusted_R_squared`: Adjusted R-squared value of the predictions.
- `Functional_Correlation`: Functional Correlation between the true and predicted response variables.
- `CDF_Distance`: Cumulative Distribution Function distance between the true and predicted response variables.
- `Area_Between_Curves`: Area between the true and predicted response variable curves.
"""
function compute_metrics(Y_test, Z_coeff, beta_matrix, beta_star, alpha_star, groups, true_predictors_train)
    Y_pred = get_predictions(Z_coeff, beta_star, alpha_star)

    # Coefficient Evaluation
    mse_beta = mse_coefficients(beta_matrix, beta_star)
    rmse_beta = rmse_coefficients(beta_matrix, beta_star)


    # Model Prediction Evaluation
    mse_Y = mse_predictions(Y_test, Y_pred)
    rmse_Y = rmse_predictions(Y_test, Y_pred)
    mae_Y = mae_predictions(Y_test, Y_pred)
    r2 = r_squared(Y_test, Y_pred)
    #  the length of true_predictors_train
    p = size(true_predictors_train)[1]
    adj_r2 = adjusted_r_squared(Y_test, Y_pred, p)

    # FDA Specific Metrics
    func_correlation = functional_correlation(Y_test, Y_pred)
    cdf_dist = cdf_distance(Y_test, Y_pred)
    area_between = area_between_curves(Y_test, Y_pred)
    ise_beta = integrated_squared_error_beta(beta_matrix, beta_star)
    iqr = calculate_robust_curve_distance_IQR(beta_matrix, beta_star, groups)

    # Return a dictionary of metrics
    Dict(
        "MSE_Coefficients" => mse_beta,
        "RMSE_Coefficients" => rmse_beta,
        "ISE_Coefficients" => ise_beta,
        "MSE_Predictions" => mse_Y,
        "RMSE_Predictions" => rmse_Y,
        "MAE_Predictions" => mae_Y,
        "R_squared" => r2,
        "Adjusted_R_squared" => adj_r2,
        "Functional_Correlation" => func_correlation,
        "CDF_Distance" => cdf_dist,
        "Area_Between_Curves" => area_between
    )
end
