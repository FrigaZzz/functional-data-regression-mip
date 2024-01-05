using JuMP
using Gurobi

function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M; intercept = false , group_limit=Inf, weights = nothing)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)
    if(weights == nothing)
        weights = ones(n)
    end
    # MIP parameters
    maxtime = 600
    out = 1
    # Create a model
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
    "TimeLimit" => 1800,  # Increased time limit to allow more thorough exploration
    "OutputFlag" => 1,  # Enable solver output for insights during solving
    "Presolve" => 0,  # Apply more presolving to simplify the model
    "Heuristics" => 0.1,  # Enable heuristics for finding good feasible solutions early
    "MIPGap" => 0.01,  # Set a smaller MIP gap for a closer optimal solution
    "Threads" => 0,  # Use all available threads for parallel computation
    "MIPFocus" => 1,  # Focus on finding feasible solutions quickly
    "NumericFocus" => 3,  # Increase numerical stability focus
    "NonConvex" => 2,  # Allow for non-convex optimization
    "OptimalityTol" => 1e-4,  # Tighten the optimality tolerance
    "IntFeasTol" => 1e-6,  # Tighten the integer feasibility tolerance
    "FeasibilityTol" => 1e-6  # Tighten the feasibility tolerance
   ))



    # Define variables
    @variable(model, beta[1:p, 1:r])
    @variable(model, beta_nonzero[1:p, 1:r], Bin)  # Binary variables indicating whether beta[j, k] is nonzero
    @variable(model, group[1:p], Bin)  # Binary variables indicating whether any coefficient in group j is nonzero


    if intercept
        @variable(model, alpha)
    else
        alpha = 0
    end

    # Assuming you have a vector 'weights' containing w_i for each observation
    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j,:] for j in 1:p))^2 for i in 1:n) +
                      lambda * sum(beta_nonzero[j, k] for j in 1:p, k in 1:r) +
                      lambda_group * sum(group[j] for j in 1:p))

    for j in 1:p
        # Constraints to link the binary variables with the beta coefficients
        #if any coefficient in a group is nonzero, that group is selected
        @constraint(model, sum(beta_nonzero[j, k] for k in 1:r) <= r * group[j])
        for k in 1:r
            # Apply Big M constraints to ensure beta values are zero if group is zero
            @constraint(model, beta[j, k]  <= BIG_M * beta_nonzero[j, k] )
            @constraint(model, beta[j, k] >= -BIG_M * beta_nonzero[j, k] ) 

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
    beta_star = beta_star .* selected
    
    return beta_star, alpha_star, selected_groups
end

