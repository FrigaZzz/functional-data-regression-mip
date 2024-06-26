using JuMP
using Gurobi
using Statistics
# Main function for mixed-integer programming functional regression
function mip_functional_regression(Y, Z, BigM,BigM_; intercept = false, group_limit=Inf, initial_beta = nothing, initial_group = nothing, initial_alpha = nothing)
    lambda = 0.01
    n, p, r = size(Z)
    group_limit = Int(min(group_limit, p))
    # Create and configure the optimization model
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
    "TimeLimit" => 1800,  # Increased time limit to allow more thorough exploration
    "OutputFlag" => 1,  # Enable solver output for insights during solving
    "Presolve" => 2,  # Apply more presolving to simplify the model
    "Heuristics" => 1,  # Enable heuristics for finding good feasible solutions early
    "MIPGap" => 0.01,  # Set a smaller MIP gap for a closer optimal solution
    "Threads" => 0,  # Use all available threads for parallel computation
    "MIPFocus" => 1,  # Focus on finding feasible solutions quickly
    "NumericFocus" => 3,  # Increase numerical stability focus
    "NonConvex" => 2,  # Allow for non-convex optimization
    "OptimalityTol" => 1e-7,  # Tighten the optimality tolerance
    "IntFeasTol" => 1e-7,  # Tighten the integer feasibility tolerance
    "FeasibilityTol" => 1e-7  # Tighten the feasibility tolerance
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
            @constraint(model, beta[j, k] <=   group[j] * 100) #BigM[j])
            @constraint(model, beta[j, k] >=  group[j] * -100)#BigM_[j])

        end
    end

    # Set initial values if provided
    if initial_beta !== nothing
        set_start_value.(beta, initial_beta)
    end
    if initial_group !== nothing
        set_start_value.(group, initial_group)
    end
    if intercept && initial_alpha !== nothing
        set_start_value(alpha, initial_alpha)
    end

    # Solve the model
    optimize!(model)

    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    group_star = JuMP.value.(group)
    alpha_star = intercept ? JuMP.value(alpha) : 0
    # # Post-process solution for binary values
    # tolerance = 1e-6
    # selected = [JuMP.value(beta[j, k]) > tolerance ? 1 : 0 for j in 1:p, k in 1:r]
    # selected_groups = [group_star[j] > tolerance ? 1 : 0 for j in 1:p]
    # beta_star = beta_star .* selected
 
 
    return beta_star, alpha_star, group_star
end

