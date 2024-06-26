using JuMP
using Gurobi
using Statistics
# Main function for mixed-integer programming functional regression
function mip_functional_regression(
            Y, Z, 
            UpperBounds, LowerBounds; 
            intercept = false, 
            group_limit=Inf, lambda = 0,
            initial_beta = nothing, initial_group = nothing, initial_alpha = nothing)
    #lambda = 10^-2  # regularization parameter, adjust as needed

    n, p, r = size(Z)

    group_limit = Int(min(group_limit, p))
    # Create and configure the optimization model
    _LPWarmStart = initial_beta !== nothing || initial_group !== nothing || initial_alpha !== nothing ? 2 : 0
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
    "TimeLimit" => 100,  # Increased time limit to allow more thorough exploration
    "OutputFlag" => 1,  # Enable solver output for insights during solving
    "Heuristics" => 1,  # Enable heuristics for finding good feasible solutions early
    "MIPGap" => 0.01,  # Set a smaller MIP gap for a closer optimal solution
    "Threads" => 16,  # Use all available threads for parallel computation
    "MIPFocus" => 1,  # Focus on finding feasible solutions quickly
    "NumericFocus" => 2,  # Increase numerical stability focus
    "NonConvex" => 2,  # Allow for non-convex optimization
    "OptimalityTol" => 1e-9,  # Tighten the optimality tolerance
    "IntFeasTol" => 1e-9,  # Tighten the integer feasibility tolerance
    "FeasibilityTol" => 1e-9,  # Tighten the feasibility tolerance,
    "LPWarmStart" => 2  # Enable warm start if initial values are provided

   ))
# Define variables
@variable(model, beta[1:p, 1:r])
@variable(model, group[1:p], Bin) # Boolean variable for group selection

alpha = intercept ? @variable(model) : 0

# Set up the objective function
if lambda>=0
    # rappresents the l2 norm of each predictor beta
    @variable(model, in_beta[1:p])
    @constraint(model, in_beta .>= 0)  # Ensure t is non-negative (absolute value)
    for i in 1:p
        @constraint(model, in_beta[i]^2 == sum(beta[i, k]^2 for k in 1:r))  # Square both sides
    end

    ## rappresents the sum of the l2 norm of each predictor betas
    @variable(model,t)
    @constraint(model, t >= 0)  # Ensure t is non-negative (absolute value)
    @constraint(model, t == sum(in_beta[j] for j in 1:p))  # Square both sides
    
    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j, :] for j in 1:p))^2 for i in 1:n) + lambda * t)

else
    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j, :] for j in 1:p))^2 for i in 1:n))
end


# Constraints to link the group selection with the beta variables
for j in 1:p
    for k in 1:r
        # Adjusted Big M constraints
        @constraint(model, beta[j, k] <= UpperBounds[j, k] * group[j])
        @constraint(model, beta[j, k] >= LowerBounds[j, k] * group[j])
    end

end

# Group limit constraint
@constraint(model, sum(group[j] for j in 1:p) == group_limit)

# Ensure exactly `group_limit` predictors (non-zero rows in `beta`)


    # Set initial values if provided
   # Set initial values if provided
    if initial_beta !== nothing
        for j in 1:p
            for k in 1:r
                JuMP.set_start_value(beta[j, k], initial_beta[j, k])
            end
        end
    end
    if initial_group !== nothing
        for j in 1:p
            JuMP.set_start_value(group[j], initial_group[j])
        end
    end
    if intercept && initial_alpha !== nothing
        JuMP.set_start_value(alpha, initial_alpha)
    end
    

    # Solve the model
    optimize!(model)

    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    group_star = round.(JuMP.value.(group))
    alpha_star = intercept ? JuMP.value(alpha) : 0
    # # Post-process solution for binary values
    beta_star = beta_star .* group_star
    
    return beta_star, alpha_star, group_star
end

