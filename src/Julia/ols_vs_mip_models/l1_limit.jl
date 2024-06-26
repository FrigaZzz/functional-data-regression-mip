using JuMP, Gurobi

function mip_functional_regression(Y, Z,  A,B; intercept=false, group_limit=Inf)
    BIG_M = 100000
    lambda = 100
    n, p, r = size(Z)
    group_limit = min(group_limit, p)
    # MIP parameters
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


    # Define variables
    @variable(model, beta[1:p, 1:r])
    @variable(model, t[1:p, 1:r] >= 0) # For L1 regularization
    @variable(model, group[1:p], Bin) # Boolean variable for group selection
    @variable(model, beta_nonzero[1:p, 1:r], Bin) # Auxiliary binary variables


    if intercept
        @variable(model, alpha)
    else
        alpha = 0
    end
    # Set up the objective function
    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j, :] for j in 1:p))^2 for i in 1:n) + lambda * sum(t[j, k] for j in 1:p, k in 1:r)
                          
    )

    # Constraints to link the group selection with the beta variables
    for j in 1:p
        for k in 1:r
            # Adjusted Big M constraints
            @constraint(model, beta[j, k] <= BIG_M * beta_nonzero[j, k])
            @constraint(model, -beta[j, k] <= BIG_M * beta_nonzero[j, k])
            @constraint(model, beta_nonzero[j, k] <= group[j])

            # Linearize the L1 norm
            @constraint(model, beta[j, k] <= t[j, k])
            @constraint(model, -beta[j, k] <= t[j, k])
        end

    end

    # Group limit constraint
    @constraint(model, sum(group) == group_limit)

    # Solve the model
    optimize!(model)


    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    alpha_star = JuMP.value.(alpha)
    group = JuMP.value.(group)

    return beta_star, alpha_star, group
end
