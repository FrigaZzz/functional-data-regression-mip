import csv
import subprocess

# Define the CSV filename
csv_filename = 'input_sim.csv'

# Read the parameters from the CSV file
def read_parameters_from_csv(csv_file):
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            yield row

# Run the Julia model execution script with the parameters
def run_julia_script(params):
    julia_command = ["julia", "model_execution.jl"] + [str(params[key]) for key in params]
    subprocess.run(julia_command, check=True)

def main():
    for params in read_parameters_from_csv(csv_filename):
        try:
            # Run the Julia script with the parameters
            run_julia_script(params)
            print(f"Executed model with parameters: {params}")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred when executing the model: {e}")

if __name__ == "__main__":
    main()
