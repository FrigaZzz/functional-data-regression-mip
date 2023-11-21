library(here)

simulation_name = "2_predictors_coef"
simulation_settings_file = "default"
# Source the generic simulator script
source(here("src", "R", "generic_simulator", "simulate_paper.R"))
source(here("src", "R", "generic_simulator", "utils","loader_utilities.R"))

set.seed(1)
# Required inputs before running the simulation!!!
inputs  <- load_simulation_settings(simulation_name, simulation_settings_file)
inputs$measurements = 10
time_domains_eval <- lapply(inputs$time_domains, function(domain) {
           seq(from = domain[[1]], to = domain[[2]], length.out = inputs$measurements)
       })
       inputs$time_domains <- time_domains_eval
inputs$noise_snr = c(100,100)

betas = create_beta_curves(beta_funcs, inputs$time_domains)
observations = 3 
# X <- simulate_functional_features(mu_funcs, cov_funcs, 1, inputs$time_domains)
U = simulate_true_predictors_Ut(mu_funcs, observations, inputs$time_domains);



basis = create_1_basis(inputs$basis_functions, inputs$time_domains[[1]], inputs$norder) 
fdPar_obj <- fda::fdPar(basis)
# B = smooth.basis(inputs$time_domains[[1]], betas[1,], fdPar_obj)$fd$coefs
# basis_values = eval.basis(inputs$time_domains[[1]], basis)
# W = smooth.basis(inputs$time_domains[[1]], U[1,1,], basis)$fd$coefs
# result <- t(W) %*% t(basis_values)
# J = inprod(basis, basis)
# zi =  t(W) %*% J
# Y_my_func= compute_Y_values(U, betas, observations = observations, predictors, time_domains = inputs$time_domains, intercept=inputs$intercept)$Y
# Y2__library = compute_Y_values(U, betas, observations = observations, predictors, time_domains = inputs$time_domains, intercept=inputs$intercept)$Y
Y = compute_Y_values_with_func(U, betas, observations = observations, predictors, time_domains = inputs$time_domains, intercept=inputs$intercept)$Y

# merge 
time_domains = inputs$time_domains
basis_functions = inputs$basis_functions
norder = inputs$norder


# 2. Create B-spline basis object
  basis_objs <- create_basis(basis_functions, time_domains, norder, predictors)
  
  # 3. Smooth beta coefficients into the B-spline basis
  result <- smooth_betas_generic(beta_funcs,basis_functions, time_domains, basis_objs)
  B <- result$beta_matrix
  basis_values <- result$basis_values
  beta_point_values <- result$beta_point_values

  # 4. Compute the J matrix (inner products for basis functions)
  J <- compute_J_matrix_generic(basis_objs, predictors, basis_functions)

  # 5. W_array: Expand X functional data into B-spline basis
  W <- compute_W_matrix_generic(U, basis_functions, time_domains, basis_objs)
#  num_observations <- dim(U)[1]
#   num_predictors <- dim(U)[2]
#   num_basis <- basis_functions
  
#   # Initialize a 3D array to store the B-spline coefficients
#   W <-  array(0, dim = c(num_observations,num_predictors,num_basis))

#   for (i in 1:num_observations) {
#     for (j in 1:num_predictors) {
#       data_vector <- U[i, j, ]
#       time <- time_domains[[j]]
#       basis_obj <- basis_objs[[j]]
#       fd_obj <- fda::smooth.basis(time, data_vector, basis_obj)$fd
#       print(fd_obj$coefs)
#       W[i,j,] <- fd_obj$coefs
#     }
#   }


  # 6. Compute Z matrix: Z = W %*% J 
  Z_matrix <- compute_Z_matrix_generic(W, J, predictors, basis_functions)