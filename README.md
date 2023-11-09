# functional-data-regression-mip
This repository presents an innovative approach to scalar-on-function linear regression analysis using a Mixed Integer Programming (MIP) solver. The method provided here addresses the challenges of functional data analysis (FDA), where the predictors are curves, images, or shapes that vary over a continuum such as time, space, or wavelength.




---
## **Automatic Environment Setups **

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
