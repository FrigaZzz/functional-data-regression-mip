
using JuMP
using Gurobi

# Main function for mixed-integer programming functional regression
function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M; intercept = false, group_limit=Inf, initial_beta = nothing, initial_group = nothing, initial_alpha = nothing)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)

    # Create and configure the optimization model
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
        "TimeLimit" => 1800,  # Increased time limit to 30 minutes
        "OutputFlag" => 1,    # Enable solver output
        "Presolve" => 0,      # Enhanced presolving
        "Heuristics" => 0,    # Enable heuristics for early feasible solutions
        "MIPGap" => 0,        # Aim for the exact optimal solution
        "Threads" => 0,       # Utilize all available threads
        "MIPFocus" => 1,      # Focus on finding feasible solutions quickly
        "NumericFocus" => 4,  # Enhanced numerical stability
        "NonConvex" => 2,     # Allow non-convex optimization
        "OptimalityTol" => 1e-6,  # Tighten optimality tolerance
        "IntFeasTol" => 1e-6,     # Tighten integer feasibility tolerance
        "FeasibilityTol" => 1e-6  # Tighten feasibility tolerance
    ))

    # Define variables for coefficients and binary indicators
    @variable(model, beta[1:p, 1:r])
    @variable(model, beta_nonzero[1:p, 1:r], Bin)  # Binary indicators for nonzero beta coefficients
    @variable(model, group[1:p], Bin)              # Binary indicators for group selection

    # Optional intercept variable
    alpha = intercept ? @variable(model) : 0

    # Objective function: Minimize squared errors with regularization terms
    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j,:] for j in 1:p))^2 for i in 1:n) +
                          lambda * sum(beta_nonzero[j, k] for j in 1:p, k in 1:r) +
                          lambda_group * sum(group[j] for j in 1:p))

    # Constraints linking beta coefficients to binary indicators
    for j in 1:p
        @constraint(model, sum(beta_nonzero[j, k] for k in 1:r) <= r * group[j])
        for k in 1:r
            # Big M method to enforce beta coefficients to zero when corresponding binary indicator is zero
            @constraint(model, beta[j, k]  <= BIG_M[j] * beta_nonzero[j, k])
            @constraint(model, beta[j, k]  >= -BIG_M[j] * beta_nonzero[j, k])
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

    # Post-process solution for binary values
    tolerance = 1e-5
    selected = [JuMP.value(beta_nonzero[j, k]) > tolerance ? 1 : 0 for j in 1:p, k in 1:r]
    selected_groups = [group_star[j] > tolerance ? 1 : 0 for j in 1:p]
    beta_star = beta_star .* selected

    return beta_star, alpha_star, selected_groups
end



"""
# Set random seed for reproducibility
using Random
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
Big_M = [100,100,100]
# now run the model

# pre compute the initial values
initial_beta = zeros(p, r).+ 10000
initial_group = zeros(p) .+ 10000
initial_alpha = 1000
beta_star, alpha_star, selected_groups = mip_functional_regression(Y, Z, 0, 0, Big_M, intercept = true, group_limit=1, initial_beta = initial_beta, initial_group = initial_group, initial_alpha = initial_alpha)

"""