using JuMP, Gurobi

function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M; intercept=false, group_limit=Inf)
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
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => maxtime,
        "OutputFlag" => out, "Presolve" => 2,
        "Heuristics" => 0, "MIPGap" => 0.000005,
        "Threads" => 1, "MIPFocus" => 0,
        "NumericFocus" => 1, "NonConvex" => 2,
        "OptimalityTol" => 0.0001, "IntFeasTol" => 1e-9,
        "FeasibilityTol" => 1e-9))

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
    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j, :] for j in 1:p))^2 for i in 1:n) + +lambda * sum(t[j, k] for j in 1:p, k in 1:r)
                           + lambda_group * sum(group[j] for j in 1:p)
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
        # Revised group selection constraint
        @constraint(model, group[j] <= sum(beta_nonzero[j, k] for k in 1:r))
    end

    # Group limit constraint
    @constraint(model, sum(group[j] for j in 1:p) <= group_limit)

    # Solve the model
    optimize!(model)


    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    alpha_star = JuMP.value.(alpha)
    group = JuMP.value.(group)

    return beta_star, alpha_star, group
end
