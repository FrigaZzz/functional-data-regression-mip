using Statistics

function integrated_squared_error_beta(beta_matrix, beta_star)
    sum((beta_matrix - beta_star) .^ 2)
end


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


function area_between_curves(Y_test, Y_pred)
    sum(abs.(Y_test - Y_pred))
end

# Measures the distance between the CDFs of the predicted and actual functions.
function cdf_distance(Y_test, Y_pred)
    # Assuming Y_test and Y_pred are sorted
    sum(abs.(cumsum(Y_test) / sum(Y_test) - cumsum(Y_pred) / sum(Y_pred)))
end

# This metric measures the correlation between the predicted and actual functions. It's particularly useful when you want to assess how well the predicted functions follow the trend of the actual functions.

function functional_correlation(Y_test, Y_pred)
    cor(Y_test, Y_pred)
end

