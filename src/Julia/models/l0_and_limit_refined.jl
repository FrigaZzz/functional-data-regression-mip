using JuMP
using Gurobi


function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M, group_limit=Inf)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)

    # MIP parameters
    maxtime = 60
    out = 1

    # Create a model with Gurobi optimizer and specific parameters
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => maxtime,
                                            "OutputFlag" => out, "Presolve" => 2,
                                            "Heuristics" => 0.5, "MIPGap" => 0.05,
                                            "Threads" => 1, "MIPFocus" => 1,
                                            "NumericFocus" => 1, "NonConvex" => 2,
                                            "OptimalityTol" => 0.01, "IntFeasTol" => 1e-6,
                                            "FeasibilityTol" => 1e-6))

    # Define variables
    @variable(model, beta[1:p, 1:r])
    @variable(model, beta_nonzero[1:p, 1:r], Bin)
    @variable(model, group[1:p], Bin)

    # Objective function:
    residuals = [sum(Z[i, :, :] .* beta) for i in 1:n]  # Compute the residuals
    prediction_error = sum((Y - residuals) .^ 2)  # Compute the error
    sparsity_term = lambda * sum(beta_nonzero)
    group_sparsity_term = lambda_group * sum(group)
    @objective(model, Min, prediction_error + sparsity_term + group_sparsity_term)

    # Constraints to link the binary variables with the beta coefficients using vectorized form
    @constraint(model, beta .<= BIG_M .* beta_nonzero)
    @constraint(model, -beta .<= BIG_M .* beta_nonzero)

    for j in 1:p
        # Link individual nonzeros to group variable
        @constraint(model, beta_nonzero[j, :] .<= group[j])
        # If any coefficient in a group is nonzero, that group is selected
        @constraint(model, sum(beta_nonzero[j, :]) <= r * group[j])
    end

    # Group-level sparsity constraint
    @constraint(model, sum(group) <= group_limit)

    # Solve the model
    optimize!(model)

    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    group_star = JuMP.value.(group)

    # Post-process the binary values to enforce 0s and 1s based on a tolerance
    tolerance = 1e-5
    selected = Int.(JuMP.value.(beta_nonzero) .> tolerance)
    selected_groups = Int.(group_star .> tolerance)
    return beta_star, selected, selected_groups
end
