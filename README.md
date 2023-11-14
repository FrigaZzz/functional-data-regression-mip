# Functional Data Regression with MIP
This repository introduces a novel approach to scalar-on-function linear regression analysis through the application of a Mixed Integer Programming (MIP) solver. This method tackles the complexities of functional data analysis (FDA), where the predictors are curves, images, or shapes that change across a continuum such as time, space, or wavelength.



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

## Julia - Manual Project Environment Setup 


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


## R - Manual Project Environment Setup 

You can either decide to install the required dependencies that are defined in the DESCRIPTION file or reproduce the whole environment, by 
installing the RENV ENVIRONMENT.

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



---

### Functional Data Simulation Process

The simulation process generates synthetic functional data by modeling both the mean function and the covariance structure for each predictor. The data is generated to have specific statistical properties that mimic real-world scenarios in functional data analysis (FDA). This simulation approach is useful for testing FDA methods and algorithms in a controlled setting where the true underlying processes are known.

#### Inputs:

1. `mu_funcs`: A list of mean functions, each representing the expected value of the data at any given time point `t`. For example, the function `sin(2 * pi * t)` implies that the mean value of the data follows a sinusoidal pattern over time.


2. `cov_funcs`: A list of covariance functions generated by `generate_covariance_function`. Each function specifies the covariance between any two time points `t` and `s`, which describes how the observations are correlated over time.Covariance functions determine how values of the function at different points in time are related to each other, with parameters controlling the overall variance (`sig2`) and the rate at which correlation between points decreases as the distance between them increases (`rho`).
3. `beta_funcs`: A list of functions representing the effect sizes over time. They are used to model the temporal dynamics of the predictors, specifying how the impact of covariates on the `Response` changes over the domain. For instance:
- A `beta_func` defined as function(t) **4 * exp(-abs(t - 1))** suggests that the effect of the predictor on the response has a peak at t = 1 and decays exponentially away from this point. They are applied to the functional predictors to simulate the dependent variable.

- A `beta_func` that returns **0 * t**, would imply that the corresponding functional predictor has no effect on the response variable.

The influence of these beta_funcs is critical when simulating the scalar response in a dataset because it defines the nature of the relationship between the predictors and the response.

4. `time`: The range of time values over which the functions are evaluated.

#### Covariance Function Generation (`generate_covariance_function`):

This function creates a covariance function based on specified parameters:
- `sig2`: The variance parameter that scales the overall level of covariance.
- `rho`: The range parameter that influences how quickly the correlation decays with time.
- `decay_type`: The type of decay function, such as "exponential" or "matern", which defines the functional form of the decay in correlation over time.

**Exponential Decay**: This type of covariance structure implies that the correlation between observations decays exponentially with the distance between time points. It is often used to model data where the correlation quickly drops off as time increases.

**Matern Decay**: This provides a more general form of covariance structure, where the rate of decay is modulated by additional parameters. It is more flexible and can model a wider range of decay behaviors.

#### Data Generation Process:

The `simulate_functional_features` function orchestrates the simulation process, using the mean and covariance functions to generate data for each functional predictor across a specified number of observations and time points.It uses `Map` to apply the simulate_data function to each set of mean and covariance functions, which generates a multi-dimensional array representing multiple functions over the time domain.

The `simulate_data` generates for each mean and covariance function, a multivariate normal dataset with mean vector `mu` (evaluated mean function at each time point) and covariance matrix `Sigma` (evaluated covariance function at each pair of time points). The result is a matrix where each row corresponds to an observation of the functional data, and each column corresponds to a time point. 

1. **Mean Function (`mu_func`)**: Determines the central tendency of the data over time. For instance, a sinusoidal mean function will result in data that oscillates in a wave-like pattern.

2. **Covariance Function (`cov_func`)**: Describes how data points are correlated over time. An exponential decay in the covariance function implies that points closer in time are more similar, with the similarity decreasing exponentially as the time gap increases.


#### Statistical Effect and Utility:

The simulated data reflects the complex interdependencies and temporal dynamics often present in real functional data. The sinusoidal means and structured covariance provide a realistic scenario where the data points are not independent but are influenced by both their position in time and their relation to other points.

By simulating data with known properties, researchers can assess the accuracy and robustness of statistical methods in uncovering these properties from observed data.

Simulating functional data is powerful for:

- Testing Statistical Methods: Researchers and analysts can assess the performance of statistical methods under controlled conditions, understanding their behavior in response to known data-generating processes.

- Model Validation: Before applying models to real-world data, which is often messy and complex, simulated data can help validate the models' assumptions and tuning parameters.

- Educational Purposes: Simulated data can help students and practitioners learn about the characteristics of functional data and the appropriate analytical techniques.

- Software Development: Developers of statistical software can use simulated datasets to test and benchmark their algorithms.

#### Utility of Beta Functions in Simulation
The use of `beta_funcs` is essential for simulating functional data in scenarios where the goal is to understand or predict a scalar outcome based on functional predictors. This is particularly useful for:

- `Developing and Testing Functional Regression Models`: By simulating data with known beta_funcs, one can assess how well a functional regression model can recover these true effects from the data.

- `Studying Variable Effects:` Researchers can simulate scenarios where the effect of a predictor varies over time or other dimensions, which is often the case in longitudinal studies or growth curve analysis.

- `Understanding Interaction Effects`: In more complex simulations, beta_funcs can be used to model interaction effects between time-varying predictors and other covariates.


--- 

