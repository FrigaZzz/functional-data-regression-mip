
using Statistics
"""
    integrated_squared_error_beta(beta_matrix, beta_star)

Calculate the integrated squared error between two beta matrices.

# Arguments
- `beta_matrix::Matrix`: A matrix of beta coefficients.
- `beta_star::Matrix`: A matrix of beta coefficients.

# Returns
- `distance::Float64`: The integrated squared error between the two beta matrices.
"""
function integrated_squared_error_beta(beta_matrix, beta_star)
    sum((beta_matrix - beta_star) .^ 2)
end


"""
    calculate_robust_curve_distance_IQR(beta_matrix, beta_star, groups)

Calculate the robust curve distance using the interquartile range (IQR) method.

# Arguments
- `beta_matrix::Matrix`: a matrix of predictor coefficients.
- `beta_star::Matrix`: a matrix of reference coefficients.
- `groups::Vector`: a vector of binary values indicating which predictors are selected.

# Returns
- `distance::Float64`: the robust curve distance calculated using the IQR method.
"""
function calculate_robust_curve_distance_IQR(beta_matrix, beta_star, groups)
    distance = 0.0
    # Loop over each predictor
    for j in 1:size(beta_matrix, 1)
        if groups[j] == 1  # Only if the predictor is selected
            predictor_coeffs = beta_matrix[j, :]
            # Calculate the interquartile range (IQR)
            q1 = quantile(predictor_coeffs, 0.25)
            q3 = quantile(predictor_coeffs, 0.75)
            iqr = q3 - q1

            # Define the non-outlier range
            non_outlier_min = q1 - 1.5 * iqr
            non_outlier_max = q3 + 1.5 * iqr

            # Compute distance only for non-outlier coefficients
            for i in 1:size(beta_matrix, 2)
                coeff = predictor_coeffs[i]
                # Include the coefficient in the distance calculation if it's not an outlier
                if coeff >= non_outlier_min && coeff <= non_outlier_max
                    distance += (coeff - beta_star[j, i])^2
                end
            end
        end
    end
    return sqrt(distance)
end


"""
    area_between_curves(Y_test, Y_pred)

Compute the area between two curves using the absolute difference between them.

# Arguments
- `Y_test::Vector{Float64}`: The true values of the curve.
- `Y_pred::Vector{Float64}`: The predicted values of the curve.

# Returns
- `Float64`: The area between the curves.
"""
function area_between_curves(Y_test, Y_pred)
    sum(abs.(Y_test - Y_pred))
end

"""
    cdf_distance(Y_test, Y_pred)

Measures the distance between the CDFs of the predicted and actual functions.

# Arguments
- `Y_test::Vector{Float64}`: Vector of actual function values.
- `Y_pred::Vector{Float64}`: Vector of predicted function values.

# Returns
- `distance::Float64`: Distance between the CDFs of the predicted and actual functions.
"""
function cdf_distance(Y_test, Y_pred)
    # Assuming Y_test and Y_pred are sorted
    sum(abs.(cumsum(Y_test, dims=1) / sum(Y_test, dims=1) - cumsum(Y_pred, dims=1) / sum(Y_pred, dims=1)))
end


"""
    functional_correlation(Y_test, Y_pred)

Compute the correlation between two sets of functional data.

# Arguments
- `Y_test::AbstractMatrix`: the true values of the functional data.
- `Y_pred::AbstractMatrix`: the predicted values of the functional data.

# Returns
- `correlation::Float64`: the correlation between `Y_test` and `Y_pred`.
"""
function functional_correlation(Y_test, Y_pred)
    cor(Y_test, Y_pred)
end

