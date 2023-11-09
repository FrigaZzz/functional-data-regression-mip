# Functional Data Regression with MIP
This repository presents an innovative approach to scalar-on-function linear regression analysis using a Mixed Integer Programming (MIP) solver. The method provided here addresses the challenges of functional data analysis (FDA), where the predictors are curves, images, or shapes that vary over a continuum such as time, space, or wavelength.



## Repository Structure

- `setup/`: Scripts and instructions for setting up the necessary environment to run the models. Follow the setup guide before proceeding with the analysis.
- `notebooks/`: A Jupyter notebook `scalar_on_function.ipynb` is provided here to demonstrate the project's functionality.
- `simulations/`: Contains R scripts that are used to generate input data for the regression models, essential for creating datasets for analysis.
- `src/`: The source code is divided into Julia and R code:
  - `Julia/models/`: Julia implementation of the MIP models using Gurobi, forming the core of the regression analysis.
  - `Julia/utils/`: Utility scripts in Julia for supporting model execution, including data manipulation and result processing.
  - `R/`: Holds R scripts for data simulation and utility functions that complement the Julia models.
    - `data_simulator/`: R scripts for simulating datasets to be used in the models.
    - `utils/`: R utility scripts such as `basis_utilities.R`, `covariance_utilities.R`, and others for tasks like plotting and model setup.
- `outputs/`: Outputs from the model runs, including:
  - `logs/`: Textual outputs like performance metrics and execution logs.
  - `plots/`: Visual outputs such as graphs and charts that illustrate the results.

## Getting Started

1. Begin by setting up the environment using the provided scripts in the `setup/` directory.
2. Explore the `notebooks/` directory to understand the practical application of the models through the `scalar_on_function.ipynb` notebook.
3. Examine the `simulations/` folder to learn about the generation of input data for the regression analysis.
4. Review the `src/` directory for the actual model implementations in Julia and supplementary R scripts.
5. Run the models and simulations, and then analyze the outputs in the `outputs/` directory to assess the models' performance.

## Outputs

The outputs of the analyses are categorized as follows:

- `logs/`: Includes CSV files detailing errors, performance metrics, and result summaries.
- `plots/`: Contains PNG images that provide visual interpretations of the results, such as predicted curves and residuals.



---
##  Environment Setups

I have created two scripts for you: one for R using `renv` and one for Julia using the standard package manager.

- Julia: install\setup_environment.jl
- R: install\setup_environment.jl

Each script will automate the setup of the project environment by installing the necessary dependencies as specified in your `renv.lock` file for R and `Project.toml` and `Manifest.toml` files for Julia.

You can also proceed with the manual configuration.

## R - Manual Project Environment Setup 


To set up the project environment and install the necessary dependencies, follow these steps:

1. **Start Julia**: Open the Julia REPL (Read-Eval-Print Loop) by accessing the Julia command-line interface. You can do this by typing `julia` in your terminal (command prompt) and pressing Enter.

2. **Activate the Project**: Once in the Julia REPL, set the active project to the one associated with the `Project.toml` and `Manifest.toml` files. Run the following Julia commands:

   ```julia
   using Pkg
   Pkg.activate(".")
   ```

   This will activate the project environment located in the current directory (denoted by `"."`), which should contain the `Project.toml` and `Manifest.toml` files.

3. **Install Dependencies**: To install all the dependencies as specified in the `Project.toml` and `Manifest.toml`, execute:

   ```julia
   Pkg.instantiate()
   ```

   This command will ensure that you are using the same versions of the packages as specified in the `Manifest.toml` file, replicating the exact project environment.

After completing these steps, you will have all the necessary packages installed, and your environment will be set up to run the project.

---


## Julia - Manual Project Environment Setup 

1. **Open R**: Navigate to the project directory and open R.

2. **Activate `renv`**: If `renv` is not already initialized for your project, do so by running:

    ```R
    renv::init()
    ```

    If `renv` is already initialized and active, you can skip the initialization step.

3. **Restore Packages**: Use the `renv::restore()` function to install the libraries from the `renv.lock` file:

    ```R
    renv::restore()
    ```

    This command will check the `renv.lock` file and install the packages as specified there.

4. **Verify Installation**: After the restoration process completes, you can verify that the packages have been installed by trying to load them in R:

    ```R
    library(fda)
    library(MASS)
    library(refund)
    ```

    If they load without error, it means the restoration process has been successful and the packages are installed.
    Remember, `renv::restore()` will install the packages into the project-specific `renv` library, not the system-wide R library. Ensure that your R session is using the correct `renv` environment by checking with `renv::paths$library()`. This ensures that when you call `library(packageName)`, it's loading the package from the project's `renv` library.
---
