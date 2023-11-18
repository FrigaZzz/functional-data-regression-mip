using JuMP
using Gurobi

function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M; intercept = false , group_limit=Inf)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)
    # MIP parameters
    maxtime = 60
    out = 1
    # Create a model
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => maxtime,
                                            "OutputFlag" => out, "Presolve" => 2,
                                            "Heuristics" => 0.5, "MIPGap" => 0.05,
                                            "Threads" => 1, "MIPFocus" => 1,
                                            "NumericFocus" => 1, "NonConvex" => 2,
                                            "OptimalityTol" => 0.01, "IntFeasTol" => 1e-6,
                                            "FeasibilityTol" => 1e-6))



    # Define variables
    @variable(model, beta[1:p, 1:r])
    @variable(model, beta_nonzero[1:p, 1:r], Bin)  # Binary variables indicating whether beta[j, k] is nonzero
    @variable(model, group[1:p], Bin)  # Binary variables indicating whether any coefficient in group j is nonzero


    if intercept
        @variable(model, alpha)
    else
        alpha = 0
    end

    @variable(model, abs_res[1:n] >= 0)
    @variable(model, delta[1:p, 2:r])
    for j in 1:p
        for k in 2:r
            @constraint(model, delta[j, k] >= beta[j, k] - beta[j, k-1])
            @constraint(model, delta[j, k] >= -(beta[j, k] - beta[j, k-1]))
        end
    end
    
    # Set up the objective function
    @objective(model, Min, sum(
        sum(abs_res[i]  for i in 1:n))  +
        lambda * sum(beta_nonzero[j, k] for j in 1:p, k in 1:r) +   # L0 norm term
        lambda * sum(sqrt(sum(beta[j, k]^2 for k in 1:r)) * group[j] for j in 1:p)+
        lambda_group * sum(delta[j, k] for j in 1:p, k in 2:r)
        )

 


    for i in 1:n
        @constraint(model, abs_res[i] >= Y[i] - alpha - sum(Z[i, :, : ] * beta[:,:]))
        @constraint(model, abs_res[i] >= -(Y[i] - alpha - sum(Z[i, :, : ] * beta[:,:])))
    end

    
    for j in 1:p
        # Constraints to link the binary variables with the beta coefficients
        #if any coefficient in a group is nonzero, that group is selected
        @constraint(model, sum(beta_nonzero[j, k] for k in 1:r) <= r * group[j])
        for k in 1:r
            # Apply Big M constraints to ensure beta values are zero if group is zero
            @constraint(model, beta[j, k] <= BIG_M * beta_nonzero[j, k])
            @constraint(model, -beta[j, k] <= BIG_M * beta_nonzero[j, k])

            # Link individual nonzeros to group variable
            @constraint(model, beta_nonzero[j, k] <= group[j])
        end
       
    end

    # Group-level sparsity constraint
    @constraint(model, sum(group[j] for j in 1:p) <= group_limit)

    # Solve the model
    optimize!(model)

    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    group_star = JuMP.value.(group)
    alpha_star = JuMP.value.(alpha)

    # Post-process the binary values to enforce 0s and 1s based on a tolerance
    tolerance = 1e-5
    selected = [JuMP.value(beta_nonzero[j, k]) > tolerance ? 1 : 0 for j in 1:p, k in 1:r]
    selected_groups = [group_star[j] > tolerance ? 1 : 0 for j in 1:p]

    return beta_star, alpha_star, selected_groups
end
