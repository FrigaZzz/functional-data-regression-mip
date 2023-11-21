using DataFrames
using CSV


"""

Parse command line arguments to extract simulation parameters.

# Arguments
- `args`: A list of command line arguments.

# Returns
- Tuple containing all parsed parameters: `model_name`, `simulation_name`, `measurements`, `observations`, `observations_test`, `basis_functions`, `error_sd`, `seed`, `λ`, `λ_group`, `M`, `override_true_predictors`.

"""
function parse_command_line_arguments(args)
    model_name = args[1]
    simulation_name = args[2]
    setting_name = args[3]  # Retrieve setting_name from the third argument
    measurements = parse(Int, args[4])
    observations = parse(Int, args[5])
    observations_test = parse(Int, args[6])
    basis_functions = parse(Int, args[7])
    # parse array that contains a boolean and a float
    noise_snr = [lowercase(args[8]) == "true", parse(Float64, args[9])]    
    seed = parse(Int, args[10])
    λ = parse(Float64, args[11])  # Assuming λ is a floating point number
    λ_group = parse(Float64, args[12])  # Assuming λ_group is a floating point number
    M = parse(Int, args[13])
    return model_name, simulation_name, setting_name, measurements, observations, observations_test, basis_functions, noise_snr, seed, λ, λ_group, M
end
"""

Ensure that a directory exists at the specified path. If the directory does not exist, it is created.

# Arguments
- `path`: Filesystem path where the directory should exist.

# Examples
```julia
julia> ensure_directory_exists("path/to/directory")
```
"""
function ensure_directory_exists(path)
    if !isdir(path)
        mkpath(path)
    end
end


"""

Define the output directory path for simulation results based on the given parameters.

# Arguments
- `simulation_name`: Name of the simulation.
- `observations`: Number of observations in the training set.
- `observations_test`: Number of observations in the test set.
- `basis_functions`: Number of basis functions used.
- `measurements`: Number of measurements.
- `lambda`: Regularization parameter lambda.
- `lambda_group`: Group-specific regularization parameter.
- `M`: Some parameter M (usage depends on context).
- `error_sd`: Standard deviation of the simulation error.
- `seed`: Seed for the random number generator.
- `project_root`: Root path for the project directory.

# Returns
- `output_dir`: The full path to the output directory.

# Examples
```julia
julia> define_output_dir("SimY", 200, 50, 5, 100, 0.01, 0.02, 10, 0.1, 12345, "path/to/project")
"path/to/project/outputs/runs/SimY/200_50_5_100_0.01_0.02_10_0.1_12345"
```
"""
function define_output_dir(model_name, simulation_name, observations, observations_test, basis_functions, measurements, lambda, lambda_group, M, error_sd, seed, project_root)
    # Define directory structure
    folder_name = "$(simulation_name)/$(model_name)/$(observations)_$(observations_test)_$(basis_functions)_$(measurements)_$(lambda)_$(lambda_group)_$(M)_$(error_sd)_$(seed)"
    output_dir = joinpath(project_root, "outputs", "runs", folder_name)
    return output_dir
end

"""
    save_model_parameters(output_dir, params)

Save the model parameters to a CSV file in the specified output directory.

# Arguments
- `output_dir`: The directory where the CSV file will be saved.
- `params`: A table or dictionary of parameters to be saved.

# Examples
```julia
julia> save_model_parameters("path/to/output", model_params)
```
"""
function save_model_parameters(output_dir, params)
    CSV.write(joinpath(output_dir, "model_params.csv"), params)
end

"""
    save_performance_evaluation(output_dir, performance_metrics)

Save performance evaluation metrics to a CSV file in the specified output directory.

# Arguments
- `output_dir`: The directory where the CSV file will be saved.
- `performance_metrics`: A dictionary or table of performance metrics to be saved.

# Examples
```julia
julia> save_performance_evaluation("path/to/output", performance_metrics)
```
"""
function save_performance_evaluation(output_dir, performance_metrics)
    performance_df = DataFrame(performance_metrics)
    CSV.write(joinpath(output_dir, "performance_evaluation.csv"), performance_df)
end

"""
    save_model_results(output_dir, results)

Save model results to a text file in the specified output directory. This includes various matrices and vectors like `beta_matrix`, `beta_star`, etc.

# Arguments
- `output_dir`: The directory where the results file will be saved.
- `results`: A dictionary containing result matrices and vectors.

# Examples
```julia
julia> save_model_results("path/to/output", simulation_results)
```
"""
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

"""
    save_model_summary(output_dir, summary)

Save a text summary of the model to a text file in the specified output directory.

# Arguments
- `output_dir`: The directory where the summary file will be saved.
- `summary`: A string containing the summary of the model.

# Examples
```julia
julia> save_model_summary("path/to/output", "Model summary text")
```
"""
function save_model_summary(output_dir, summary)
    summary_file_path = joinpath(output_dir, "model_summary.txt")
    open(summary_file_path, "w") do io
        write(io, summary)
    end
end


"""
    print_matrix_or_tuple(io, data, name)

Prints the contents of `data` to the given `io` stream. If `data` is a tuple, it iterates through the tuple and prints each matrix element prefixed with `name` and the element index. If `data` is a single matrix, it prints the matrix prefixed with `name`.

# Arguments
- `io`: The IO stream to which the data will be printed.
- `data`: The tuple or matrix to be printed. If a tuple, each matrix element is printed separately.
- `name`: A descriptive name used as a prefix when printing the elements of `data`.
"""
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

"""

Prints a matrix with each row on a new line.

# Arguments
- `io`: An IO stream to write the output to.
- `matrix`: The matrix to be printed.
- `matrix_name`: The name of the matrix to be printed.

"""
function print_matrix(io, matrix, matrix_name)
    write(io, "$matrix_name:\n")
    for row in eachrow(matrix)
        write(io, join(row, ", "), "\n")
    end
    write(io, "\n")
end


"""
Function to read model results

Reads the results of a model from a file located at `output_dir` and returns them as a dictionary of matrices.
Each matrix is identified by a string key, which is the name of the matrix in the file.
The file must have the following format:
- Each matrix is preceded by a line containing the matrix name followed by a colon.
- Each row of the matrix is on a separate line.
- The elements of each row are separated by commas and optionally spaces.

# Arguments
- `output_dir`: The directory containing the model results file.
"""
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

