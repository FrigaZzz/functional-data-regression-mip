
import csv
import subprocess
import os

# Define the CSV filename
csv_filename = 'input_sim.csv'

# Read the parameters from the CSV file


def read_parameters_from_csv(csv_file):
    with open(csv_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            yield row


def generate_sh_script(params, time="0-01:10", mem="80", cpus_per_task="20"):
    # Construct the command to run 'model_execution.jl' with the parameters from the CSV
    julia_command = f"srun julia model_execution.jl {params['model_name']} {params['simulation_name']} {params['observations']} {params['observations_test']} {params['measurements']} {params['basis_functions']} {params['error_sd']} {params['seed']} {params['λ']} {params['λ_group']} {params['M']}"

    # SBATCH script
    return f"""#!/bin/bash
            #SBATCH --job-name={params['simulation_name']}
            #SBATCH --time={time}
            #SBATCH --mem={mem}
            #SBATCH --cpus-per-task={cpus_per_task}
            #SBATCH --output=output/runs/{params['simulation_name']}/log/Z_%A_%a.out.log
            #SBATCH --error=output/runs/{params['simulation_name']}/log/Z_%A_%a.out.log

            # Load Julia module if needed
            module load gurobi/9.5.0
            module load julia/1.9.1
            module load r/4.3.1

            export seed=$SLURM_ARRAY_TASK_ID

            # Run the model
            {julia_command} > output/runs/{params['simulation_name']}/log/${{SLURM_ARRAY_JOB_ID}}_{params['observations']}_{params          ['observations_test']}_{params['basis_functions']}_{params['measurements']}_{params['λ']}_{params['λ_group']}_{params['M']}_{params         ['error_sd']}_{params['seed']}.log
            """


def main():
    for params in read_parameters_from_csv(csv_filename):
        sh_script = generate_sh_script(params)
        script_filename = f"{params['simulation_name']}.sh"

        # Write the shell script to a file
        with open(script_filename, 'w') as script_file:
            script_file.write(sh_script)

        try:
            # Make the script executable
            subprocess.run(['chmod', '+x', script_filename], check=True)

            # Submit the script to SBATCH
            sbatch_command = f"sbatch --array=1-10%1 {script_filename}"
            subprocess.run(sbatch_command, check=True, shell=True)
            print(f"Submitted job array with parameters: {params}")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred when submitting the job array: {e}")

        # Delete the script file
        os.remove(script_filename)


if __name__ == "__main__":
    main()
