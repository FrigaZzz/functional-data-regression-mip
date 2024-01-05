# import csv
# import subprocess

# # Define the CSV filename
# csv_filename = 'input_sim.csv'

# # Read the parameters from the CSV file
# def read_parameters_from_csv(csv_file):
#     with open(csv_file, newline='') as csvfile:
#         reader = csv.DictReader(csvfile)
#         for row in reader:
#             yield row

# # Run the Julia model execution script with the parameters
# def run_julia_script(params):
#     julia_command = ["julia", "model_execution.jl"] + [str(params[key]) for key in params]
#     subprocess.run(julia_command, check=True)

# def main():
#     for params in read_parameters_from_csv(csv_filename):
#         try:
#             # Run the Julia script with the parameters
#             run_julia_script(params)
#             print(f"Executed model with parameters: {params}")
#         except subprocess.CalledProcessError as e:
#             print(f"An error occurred when executing the model: {e}")

# if __name__ == "__main__":
#     main()

import csv
import subprocess
import concurrent.futures

# Define the CSV filename
csv_filename = 'input_sim.csv'

# Number of rows per thread
rows_per_thread = 10

# Read the parameters from the CSV file
def read_parameters_from_csv(csv_file):
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        batch = []
        for row in reader:
            batch.append(row)
            if len(batch) >= rows_per_thread:
                yield batch
                batch = []
        if batch:  # Handle any remaining rows
            yield batch

# Run the Julia model execution script with the parameters
def run_julia_script(params):
    julia_command = ["julia", "model_execution.jl"] + [str(params[key]) for key in params]
    subprocess.run(julia_command, check=True)

# Function to process a batch of rows
def process_batch(batch):
    for params in batch:
        try:
            # Run the Julia script with the parameters
            run_julia_script(params)
            print(f"Executed model with parameters: {params}")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred when executing the model: {e}")

def main():
    # Using ThreadPoolExecutor to parallelize the task
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for batch in read_parameters_from_csv(csv_filename):
            executor.submit(process_batch, batch)

if __name__ == "__main__":
    main()
