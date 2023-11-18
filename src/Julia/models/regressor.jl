
using JuMP
using Gurobi

function mip_functional_regression(Y, Z, lambda, lambda_group, BIG_M; intercept=false, group_limit=Inf)
    n, p, r = size(Z)
    group_limit = min(group_limit, p)
    # MIP parameters
    maxtime = 60
    out = 1

    # outliers
    frcont = 0.35
    k_n = floor(n * frcont)

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


    
    # Set up the objective function
    @objective(model, Min, sum((Y[i]  - sum(Z[i, :, : ] * beta[:,:]) )^2 for i in 1:n))



    # Solve the model
    optimize!(model)

    # Retrieve the solution
    beta_star = JuMP.value.(beta)
    


    return beta_star, 0,0
end
