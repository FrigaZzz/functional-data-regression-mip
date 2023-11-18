using JuMP, Gurobi

function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M; intercept = false , group_limit=Inf)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)
    # MIP parameters
    MIPpresolve = 2
    MIPheur = 0.5
    MIPfocus = 1
    maxtime = 60
    MIPg = 0.05
    threads = 1
    out = 1

    # Create a model
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => maxtime, "Heuristics" => MIPheur,
        "MIPGap" => MIPg, "Threads" => threads, "OutputFlag" => out, "Presolve" => MIPpresolve,
        "MIPfocus" => MIPfocus, "NumericFocus" => 1, "NonConvex" => 2, "OptimalityTol" => 0.01,
        "IntFeasTol" => 1e-9,"FeasibilityTol" =>1e-6))

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
    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, : ]' * beta[j,:] for j in 1:p))^2 for i in 1:n) +

        + lambda * sum(t[j, k] for j in 1:p, k in 1:r)
        + lambda_group * sum(group[j] for j in 1:p)
    )

    # Constraints to link the group selection with the beta variables
    for j in 1:p
        # Enforce that if any beta for a group is nonzero, the group must be selected
        @constraint(model, sum(beta_nonzero[j, k] for k in 1:r) >= group[j])
        for k in 1:r
            # Apply Big M constraints to ensure beta values are zero if group is zero
            @constraint(model, beta[j, k] <= BIG_M * beta_nonzero[j, k])
            @constraint(model, -beta[j, k] <= BIG_M * beta_nonzero[j, k])

            # Enforce that if a group is selected, at least one beta must be nonzero
            @constraint(model, beta_nonzero[j, k] <= group[j])
            
            # Linearize the L1 norm
            @constraint(model, beta[j, k] <= t[j, k])
            @constraint(model, -beta[j, k] <= t[j, k])
        end
    end

    # Constraint to enforce sparsity at the group level (optional, can be adjusted)
    # This constraint ensures that the sum of all group variables is less than or equal to a certain limit
    @constraint(model, sum(group[j] for j in 1:p) <= group_limit) # where group_limit is a parameter you can define


    # Solve the model
    optimize!(model)

    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    z_star = JuMP.value.(Z)
    alpha_star = JuMP.value.(alpha)
    group = JuMP.value.(group)

    # Post-process the binary values to enforce 0s and 1s based on a tolerance
    tolerance = 1e-5
    selected = [value > tolerance ? 1 : 0 for value in z_star]
    return beta_star, alpha_star, group
end
