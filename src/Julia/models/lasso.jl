using Gurobi, JuMP

function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M; intercept = false , group_limit=Inf)
    n, p, r = size(Z)

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
    "MIPfocus" => MIPfocus, "NumericFocus"=>1, "NonConvex"=>2, "OptimalityTol" => 0.01))

    # Define variables
    @variable(model, beta[1:p, 1:r])
    @variable(model, t[1:p, 1:r] >= 0) # For L1 regularization
    @variable(model, alpha[1:p])
    @variable(model, group[1:p], Bin) # Boolean variable for 

    # Set up the objective function
    @objective(model, Min, 
        sum((Y[i] - sum(Z[i, :, : ] * beta[:,:]))^2 for i in 1:n) 
        + lambda * sum(t[j, k] for j in 1:p, k in 1:r) +
        + lambda_group * sum(group[j] for j in 1:p)
    )

     for j in 1:p
        for k in 1:r
            
            # Apply Big M constraints to ensure beta values are zero if group is zero
            @constraint(model, beta[j, k] <= M * group[j])
            @constraint(model, beta[j, k] >= -M * group[j])

            # Linearize the L1 norm
            @constraint(model, beta[j, k] <= t[j, k])
            @constraint(model, -beta[j, k] <= t[j, k])

        end
    end

    # Solve the model
    optimize!(model)

    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    alpha_star = JuMP.value.(alpha)
    selected_raw = JuMP.value.(group)

    # Post-process the group values to enforce 0s and 1s based on a tolerance
    tolerance = 1e-5
    selected = [value > tolerance ? 1 : 0 for value in selected_raw]
    
    return beta_star, alpha_star, selected_raw
end
