using JuMP
using Gurobi

function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_MU, BIG_ML; intercept = false , group_limit=Inf, initial_beta = nothing, initial_group = nothing, initial_alpha = nothing)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)
    # MIP parameters
    maxtime = 600
    out = 1
    # Create a model
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
    "OptimalityTol" => 1e-5,  # Tighten the optimality tolerance
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
    @objective(model, Min,1/n *sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j,:] for j in 1:p))^2 for i in 1:n)  +
    lambda * sum(beta_nonzero[j, k] for j in 1:p, k in 1:r) +
    lambda_group * sum(group[j] for j in 1:p))

    
    for j in 1:p
        # Constraints to link the binary variables with the beta coefficients
        #if any coefficient in a group is nonzero, that group is selected
        @constraint(model, sum(beta_nonzero[j, k] for k in 1:r) <= r * group[j])
        for k in 1:r
            # Apply Big M constraints to ensure beta values are zero if group is zero
            @constraint(model, beta[j, k]  <= BIG_MU[j,k] * beta_nonzero[j, k] )
            @constraint(model, beta[j, k]  >= BIG_ML[j,k] * beta_nonzero[j, k] ) 

        end
       
    end

    # Group-level sparsity constraint
    @constraint(model, sum(group[j] for j in 1:p) <= group_limit)
    
    # Set initial values if provided
    initial_beta !== nothing && set_start_value.(beta, initial_beta)
    initial_group !== nothing && set_start_value.(group, initial_group)
    intercept && initial_alpha !== nothing && set_start_value(alpha, initial_alpha)

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

"""


# Set random seed for reproducibility
Random.seed!(1234)

# test the model on a simulated dataset
n = 100
p = 3
r = 4
Z = zeros(n, p, r)

# only the first group is nonzero and follows a normal distribution with mean 1 and sd 1
Z[:, 1, 1] = randn(n) .+ 1
Z[:, 1, 2:end] = randn(n, r-1)

Z[:, 2, 1] = randn(n) .+ 1
Z[:, 2, 2:end] = randn(n, r-1)

# Define a specific weight matrix B
B = randn(p, r)
# Set the second and third group to zero
B[2:end, :] = B[2:end, :]  .* 0 

# Generate Y using the weight matrix B
Y = zeros(n)
for i in 1:n
    Y[i] = sum(Z[i, :, :] .* B)
end

# Add noise to Y
# Y = Y .+ randn(n)

# now run the model
beta_star, alpha_star, selected_groups = mip_functional_regression(Y, Z, 0, 0, 1000000, intercept = true, group_limit=1)

"""