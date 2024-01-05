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

### Julia - Manual Project Environment Setup 


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


### R - Manual Project Environment Setup 

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




## Functional Data Simulation Pipeline

The main execution pipeline involves several steps, from generating functional predictors (both true and observed) to computing response variables (Y) and incorporating errors and corrective terms.
Additional steps include creating B-spline basis objects, smoothing beta coefficients, and computing matrices essential for functional data analysis (e.g., \( J \), \( W \), \( Z \) matrices).

## Generation of True and Observed Data

### True Data Generation (`simulate_true_predictors_Ut`)

- **Purpose**: Generates true functional predictors (\( U \)) based on provided mean functions (\( mu\_funcs \)), number of observations, and time domains.
- **Mean Functions (\( mu\_funcs \))**: These are defined to allow for variation across observations. Each mean function can incorporate `random coefficients` drawn from distributions like normal or uniform or through `a matern process`. This leads to each observation having a unique realization of the mean function. 

### Observed Data Generation (`simulate_observations_Xt`)

- **Purpose**: Applies amplitude normalization to the true data \( U \) and adds random error terms to generate observed functional predictors (\( X \)).
- **Process**:
  - For each observation and predictor, a range-based standard deviation of the error term is calculated.
  - Random error terms are generated and added to \( U \) to form \( X \), introducing variability and mimicking observational noise.

## Response Variable (Y) Generation and Error Addition

### Computing Y Values (`compute_Y_values`)

- **Purpose**: Computes the response variable \( Y \) based on the observed functional predictors \( X \), beta coefficients, and other parameters like the intercept.
- **Method**:
  - Functional data objects are created for both \( X \) and the beta coefficients.
  - These objects are multiplied and integrated over the time domain to compute \( Y \).

### Amplitude Normalization and Error Addition

- **Amplitude Normalization (`compute_amplitude_norm`)**: This function computes a normalization factor based on the range of \( Y \) and adds normally distributed random errors, introducing an amplitude-related variability to \( Y \).
- **Signal-to-Noise Ratio (SNR) Adjustment (`add_snr_noise`)**: Optionally, noise can be added to \( Y \) to achieve a desired SNR level, providing control over the signal quality.


## Functional Data Simulation Strategies

We will now outlinetwo the distinct strategies for simulating functional data. Each strategy employs a different approach to the generation of mean functions, impacting how the functional data varies across observations.

### Strategy 1: Consistent Mean Across Observations

This strategy involves simulating functional data where all observations share the same mean function for each feature. The variability among observations is introduced through a `covariance function`, leading to different realizations around a shared mean for a `specific predictor`.

This strategy is appropriate when the underlying mean behavior of the feature is consistent across all observations, but individual variations are expected due to random effects.

#### Inputs:


1. `mu_funcs`
   -  A list of mean functions, each representing the expected value of the data at any given time point `t`. 
    - They are deterministic and consistent across all observations, defining a common mean trend or pattern for each feature.
    - They are defined over a `specific time domain`
```r
  mu_funcs <- list(
    function(t) sin(2 * pi * t),
    function(t) 4 * t * (1 + cos(2 * pi * t)),
    function(t) 4 * t^2
  )
```
2. `time_domains`: The range of time domains over which each function is evaluated.
```r
  time_domains = list(
  list(0, 1),
  list(0, pi / 3),
  list(-1, 1)
  )
```
3. `cov_funcs`: 
   - A list of covariance functions generated by `generate_covariance_function`. 
   - Each function specifies the covariance between any two time points `t` and `s`, which describes how the observations are correlated over time.
   - Covariance functions determine how values of the function at different points in time are related to each other, with parameters controlling the overall variance (`sig2`) and the rate at which correlation between points decreases as the distance between them increases (`rho`).
```r
     cov_funcs <- list(
    list(sig2 = 0.5, rho = 1, decay_type = "matern"),
    list(sig2 = 0.5, rho = 1, decay_type =  "matern"),
    list(sig2 = 0.5, rho = 1, decay_type =  "matern")
  )
```

4. `beta_funcs`: A list of functions representing the effect sizes over time. They are used to model the temporal dynamics of the predictors, specifying how the impact of covariates on the `Response` changes over the domain. For instance:
```r
beta_funcs <- list(
  function(t) {
    4 * exp(-abs(t - 1))
  },
  function(t) {
    sin(2 * t)
  },
  function(t) {
    rep(0, length(t))
  }
)
```
- A `beta_func` defined as function(t) **4 * exp(-abs(t - 1))** suggests that the effect of the predictor on the response has a peak at t = 1 and decays exponentially away from this point. They are applied to the functional predictors to simulate the dependent variable.

- A `beta_func` that returns **rep(0, length(t))**, would imply that the corresponding functional predictor has no effect on the response variable.

The influence of these beta_funcs is critical when simulating the scalar response in a dataset because it defines the nature of the relationship between the predictors and the response.

#### Covariance Function Generation (`generate_covariance_function`):

This function creates a covariance function based on specified parameters:
- `sig2`: The variance parameter that scales the overall level of covariance.
- `rho`: The range parameter that influences how quickly the correlation decays with time.
- `decay_type`: The type of decay function, such as "exponential" or "matern", which defines the functional form of the decay in correlation over time.

**Exponential Decay**: This type of covariance structure implies that the correlation between observations decays exponentially with the distance between time points. It is often used to model data where the correlation quickly drops off as time increases.

**Matern Decay**: This provides a more general form of covariance structure, where the rate of decay is modulated by additional parameters. It is more flexible and can model a wider range of decay behaviors.

#### Data Generation Process:

The `simulate_functional_features` function orchestrates the simulation process, using the mean and covariance functions to generate data for each functional predictor across a specified number of observations and time points. It generates data for each feature by drawing from a multivariate normal distribution with the mean provided by mu_funcs and the covariance matrix provided by the covariance function (`e.g. Matern Process`).

The `simulate_data` generates for each mean and covariance function, a multivariate normal dataset with mean vector `mu` (evaluated mean function at each time point) and covariance matrix `Sigma` (evaluated covariance function at each pair of time points). 
This results in each feature having its own characteristic mean behavior (as defined by mu_funcs) and variability and correlation structure (as defined by the Matérn covariance function).

1. **Mean Function (`mu_func`)**: Determines the central tendency of the data over time. For instance, a sinusoidal mean function will result in data that oscillates in a wave-like pattern.

2. **Covariance Function (`cov_func`)**: Describes how data points are correlated over time. An exponential decay in the covariance function implies that points closer in time are more similar, with the similarity decreasing exponentially as the time gap increases.


Using `MASS::mvrnorm` to simulate the data can address the issue of introducing variability across observations for each predictor. 

1. **Covariance Structure**: By using `mvrnorm`, one can specify a covariance matrix that describes the relationships between the predictors. This way, one can simulate data such that the values for each predictor are not only random but also correlated in a specific way. This is essential if the predictors are expected to have some degree of correlation (e.g., measurements that are physiologically linked). In our case, it is important that each feature is considered indipendent by any other.

1. **Sampling**: Instead of setting the coefficients at the beginning and keeping them constant across all observations, `mvrnorm` allows one to sample a new set of coefficients for each observation based on the multivariate normal distribution. This means that for each observation, while the mean of the coefficients remains the same (as defined by the model), their actual values will differ due to the randomness introduced by the sampling process.

3. **Variability Across Observations**: It leads to variability across observations for each predictor. This models the real-world scenario where different instances or subjects have similar but not identical characteristics.

### Strategy 2: Variable Mean Functions Across Observations

In contrast to Strategy 1, this approach allows the mean function itself to vary for each observation. This is achieved by incorporating parameters into the mean functions that are specific to each observation.

This strategy is suitable for scenarios where the mean behavior of the feature is expected to vary significantly among observations, capturing more complex and individual-specific patterns.
The main implementation is inside the function `simulate_true_predictors_Ut`

#### Mean Functions with Parameters
The mean functions are defined with parameters that vary for each observation. For example:

```r
mu_funcs <- list(
  function(t, args) { cos(2 * pi * (t - args$a1)) + args$a2 },
  function(t, args) { args$b1 * sin(pi * t) + args$b2 },
  function(t, args) { args$c1 * t^3 + args$c2 * t^2 + args$c3 }
  # Additional functions as needed
)
```

The coefficients (`a1`, `a2`, `b1`, `b2`, etc.) are generated for each predictor and each observation, typically drawn from a normal or uniform distribution.
For each observation and each predictor, the function applies the mean function with the specific coefficients, resulting in `a unique mean function for each observation`.

Same definitions apply for the time domains, beta_funcs objects.



