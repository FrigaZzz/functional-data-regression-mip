using DataFrames
using CSV


# Function to retrieve and parse command line arguments
function parse_command_line_args(args)
    model_name = args[1]
    simulation_name = args[2]
    measurements = parse(Int, args[3])
    observations = parse(Int, args[4])
    observations_test = parse(Int, args[5])
    basis_functions = parse(Int, args[6])
    error_sd = parse(Float64, args[7])
    seed = parse(Int, args[8])
    λ = parse(Float64, args[9])  # Assuming λ is a floating point number
    λ_group = parse(Float64, args[10])  # Assuming λ_group is a floating point number
    M = parse(Int, args[11])
    if length(args) > 11
        override_true_predictors = parse(Int, args[12])
    else
        override_true_predictors = nothing
    end
    return model_name, simulation_name, measurements, observations, observations_test, basis_functions, error_sd, seed, λ, λ_group, M, override_true_predictors 
end


# Function to ensure a directory exists
function ensure_directory_exists(path)
    if !isdir(path)
        mkpath(path)
    end
end

function define_output_dir(simulation_name, observations, observations_test, basis_functions, measurements, lambda, lambda_group, M, error_sd, seed, project_root)
    # Define directory structure
    folder_name = "$(simulation_name)/$(observations)_$(observations_test)_$(basis_functions)_$(measurements)_$(lambda)_$(lambda_group)_$(M)_$(error_sd)_$(seed)"
    output_dir = joinpath(project_root, "outputs", "runs", folder_name)
    return output_dir
end

# Function to save model parameters
function save_model_parameters(output_dir, params)
    params_df = DataFrame(params)
    CSV.write(joinpath(output_dir, "model_params.csv"), params_df)
end

# Function to save performance evaluation
function save_performance_evaluation(output_dir, performance_metrics)
    performance_df = DataFrame(performance_metrics)
    CSV.write(joinpath(output_dir, "performance_evaluation.csv"), performance_df)
end

# Function to save model results
function save_model_results(output_dir, results)
    results_file_path = joinpath(output_dir, "model_results.txt")
    open(results_file_path, "w") do io
        # access the betastar key in results dictionary
        print_matrix_or_tuple(io, results["beta_matrix"][1], "Beta Matrix")
        print_matrix_or_tuple(io, results["beta_star"][1], "Beta Star")
        print_matrix(io, results["true_predictors"][1], "True Predictors")
        print_matrix(io, results["groups"][1], "Evaluated Predictors")
        print_matrix(io, results["True Intercept"][1], "True Intercept")
        print_matrix(io, results["Evaluated Intercept"][1], "Evaluated Intercept")
    end
end

# Function to save model summary
function save_model_summary(output_dir, summary)
    summary_file_path = joinpath(output_dir, "model_summary.txt")
    open(summary_file_path, "w") do io
        write(io, summary)
    end
end


# Updated Helper function to print a matrix or a tuple of matrices
function print_matrix_or_tuple(io, data, name)
    if data isa Tuple
        for (index, element) in enumerate(data)
            if element isa AbstractMatrix
                print_matrix(io, element, "$(name) - Element $index")
            end
        end
    elseif data isa AbstractMatrix
        print_matrix(io, data, name)
    end
end

# Function to print a matrix with each row on a new line
function print_matrix(io, matrix, matrix_name)
    write(io, "$matrix_name:\n")
    for row in eachrow(matrix)
        write(io, join(row, ", "), "\n")
    end
    write(io, "\n")
end


# Function to read model results
function read_model_results(output_dir)
    results_file_path = output_dir
    results = Dict()
    current_key = ""
    open(results_file_path, "r") do io
        for line in eachline(io)
            if endswith(line, ":")
                # This is a matrix name
                current_key = strip(chop(line, tail=1)) # Remove the trailing colon
                results[current_key] = []
            elseif !isempty(line)
                # This is a row of a matrix
                row = parse.(Float64, split(line, ", "))
                push!(results[current_key], row)
            end
        end
    end
    # Convert lists of rows to matrices
    for (key, value) in results
        results[key] = hcat(value...)'
    end
    return results
end

