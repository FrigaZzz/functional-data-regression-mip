using JuMP
using Gurobi
using Statistics
# Main function for mixed-integer programming functional regression
function mip_functional_regression(Y, Z, BigM,BigM_; intercept = false, group_limit=Inf, initial_beta = nothing, initial_group = nothing, initial_alpha = nothing, lambda = 0)
    
    n, p, r = size(Z)
    group_limit = Int(min(group_limit, p))
    # Create and configure the optimization model
    _LPWarmStart = initial_beta !== nothing || initial_group !== nothing || initial_alpha !== nothing ? 2 : 0
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
    "TimeLimit" => 20,  # Increased time limit to allow more thorough exploration
    "OutputFlag" => 1,  # Enable solver output for insights during solving
    "Heuristics" => 1,  # Enable heuristics for finding good feasible solutions early
    "MIPGap" => 0.01,  # Set a smaller MIP gap for a closer optimal solution
    "Threads" => 0,  # Use all available threads for parallel computation
    "MIPFocus" => 1,  # Focus on finding feasible solutions quickly
    "NumericFocus" => 2,  # Increase numerical stability focus
    "NonConvex" => 2,  # Allow for non-convex optimization
    "OptimalityTol" => 1e-9,  # Tighten the optimality tolerance
    "IntFeasTol" => 1e-9,  # Tighten the integer feasibility tolerance
    "FeasibilityTol" => 1e-9,  # Tighten the feasibility tolerance,
    "LPWarmStart" => 2  # Enable warm start if initial values are provided

   ))
    # model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
    # "TimeLimit" => 1800, "OutputFlag" => 1, "Threads" => 0))


    @variable(model, beta[1:p, 1:r])
    @variable(model, group[1:p], Bin)
    alpha = intercept ? @variable(model) : 0

    @objective(model, Min, mean(sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j, :] for j in 1:p))^2 for i in 1:n)))

    # # Constraint to limit the number of groups
    @constraint(model, sum(group) == group_limit)

    # Constraints to link group and beta
    for j in 1:p
        for k in 1:r
            @constraint(model, beta[j, k] <=   group[j] * BigM[j,k])
            @constraint(model, beta[j, k] >=  group[j] * BigM_[j,k])

        end
    end

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

