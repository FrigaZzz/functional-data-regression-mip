using JuMP
using Gurobi
using Statistics
# Main function for mixed-integer programming functional regression
# Main function for mixed-integer programming functional regression
function mip_functional_regression(Y, Z, BigM, BigM_; intercept=false, group_limit=Inf, initial_beta=nothing, initial_group=nothing, initial_alpha=nothing, lambda=0.01)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)

#     # Create and configure the optimization model
#     model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
#     "TimeLimit" => 1800,  # Increased time limit to allow more thorough exploration
#     "OutputFlag" => 1,  # Enable solver output for insights during solving
#     "Presolve" => 2,  # Apply more presolving to simplify the model
#     "Heuristics" => 1,  # Enable heuristics for finding good feasible solutions early
#     "MIPGap" => 0.01,  # Set a smaller MIP gap for a closer optimal solution
#     "Threads" => 0,  # Use all available threads for parallel computation
#     "MIPFocus" => 1,  # Focus on finding feasible solutions quickly
#     "NumericFocus" => 3,  # Increase numerical stability focus
#     "NonConvex" => 2,  # Allow for non-convex optimization
#     "OptimalityTol" => 1e-5,  # Tighten the optimality tolerance
#     "IntFeasTol" => 1e-6,  # Tighten the integer feasibility tolerance
#     "FeasibilityTol" => 1e-6  # Tighten the feasibility tolerance
#    ))
    model = Model(optimizer_with_attributes(Gurobi.Optimizer,
        "TimeLimit" => 1800, "OutputFlag" => 1, "Threads" => 0))


    @variable(model, beta[1:p, 1:r])
    @variable(model, group[1:p], Bin)

    @variable(model, tbeta[1:p, 1:r])

    alpha = intercept ? @variable(model) : 0

    @objective(model, Min, sum((Y[i] - alpha - sum(Z[i, j, :]' * beta[j, :] for j in 1:p))^2 for i in 1:n) + lambda * sum(tbeta[j, k] for j in 1:p, k in 1:r) / n)



    # # Constraint to limit the number of groups
    @constraint(model, sum(group) <= group_limit)
    # @constraint(model, sum(group) >= group_limit)


    # Constraints to link group and beta
    for j in 1:p
        for k in 1:r
            @constraint(model, beta[j, k] <= group[j] * BigM[j])
            @constraint(model, beta[j, k] >= group[j] * BigM_[j])

            # Linearize the L1 norm
            @constraint(model, beta[j, k] <= tbeta[j, k])
            @constraint(model, -beta[j, k] <= tbeta[j, k])
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
    # tolerance = 1e-5
    # selected = [JuMP.value(beta[j, k]) > tolerance ? 1 : 0 for j in 1:p, k in 1:r]
    # selected_groups = [group_star[j] > tolerance ? 1 : 0 for j in 1:p]
    # beta_star = beta_star .* selected

    return beta_star, alpha_star, group_star
end


# # @__DIR__ is the directory of the current file
# # We need to go up to the parent directory to find the project root
# project_root = dirname(dirname(dirname(@__DIR__)))


# include(joinpath(project_root, "setup", "init_env.jl"))
# set_R_lib_path(project_root)


# include(joinpath(project_root, "src", "Julia", "utils", "simulation.jl"))


# simulation_name = "robust"
# simulation_settings_file = "default"

# params_train = (
#     observations = 600,
#     measurements = 300,
#     basis_functions = 12,
#     seed = 1
# )

# params_test = (
#     observations = 200,
#     measurements = 300,
#     basis_functions = 12,
#     seed = 22
# )

# # Note: Use ... to unpack NamedTuple into keyword arguments
# output = load_simulation_robust(simulation_name, simulation_settings_file, project_root; params_train...)

# output_test = load_simulation_robust(simulation_name, simulation_settings_file, project_root; params_test...)


# # Grab the outputs from the R script

# predictors = Int(output[:predictors])
# true_predictors = output[:true_predictors]
# intercept = output[:intercept]
# observations = Int(output[:observations])

# # betas and basis
# beta_matrix  = output[:B]
# basis_objs   = output[:basis_objs]
# basis_values = output[:basis_values]
# time_domains = output[:time_domains]

# # matrixes 
# U = output[:U]
# X = output[:X]
# Y = output[:Y]
# Z = output[:Z]
# J = output[:J]
# W = output[:W]


# X_test = output_test[:X]
# Y_test = output_test[:Y]
# Z_test = output_test[:Z]
# J_test = output_test[:J]
# W_test = output_test[:W]




# beta_matrix_max_values = maximum(beta_matrix, dims = 2)
# beta_matrix_min_values = minimum(beta_matrix, dims = 2)

# BigM = ones(predictors) .*  2000.0
# BigM_ =  ones(predictors) .* -2000.0

# group_limit = sum(true_predictors)
# print(group_limit)
# intercept = output[:intercept] != 0
# beta_star, alpha_star, groups = mip_functional_regression(Y, Z, BigM,BigM_; intercept = output[:intercept] != 0, group_limit= to_predict)


