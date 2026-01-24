# ==============================================================================
# 6_2_beta_estimation.R
# Section 6.2
# Estimates Beta for a SINGLE hyperparameter pair.
# ==============================================================================
# Installation of FSHybridPLS package
# This package is generated using litr, a literate programming tool developed by Jacob Bien and Patrick Vossler (https://jacobbien.github.io/litr-project/). I used litr to build the package and you don't have to install litr just to use the package. To install the package:
#
# devtools::install("your_downloaded_directory/FSHybridPLS")

library(fda)
library(dplyr)
library(tidyr)
library(FSHybridPLS)

# ------------------------------------------------------------------------------
# 0. SETUP & ARGUMENTS
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 0. SETUP & ARGUMENTS
# ------------------------------------------------------------------------------
# Default User Configuration
REPO_DIR <- "."
SAVE_DIR_REL <- "./results/beta_estimation"
REP_ID_DEFAULT <- 1
N_SAMPLE_DEFAULT <- 1000
L1_DEFAULT <- 0.1
L2_DEFAULT <- 0.01

# Create 'results' directory if it doesn't exist
if (!dir.exists(file.path(REPO_DIR, "results"))) dir.create(file.path(REPO_DIR, "results"))
save_dir <- file.path(REPO_DIR, SAVE_DIR_REL)
if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 4) {
  # 1. Parse Replication ID
  rep_id <- as.integer(args[1])
  # 2. Parse Sample Size
  n_sample_input <- as.integer(args[2])
  # 3. Parse Specific Lambda Values
  l1_input <- as.numeric(args[3])
  l2_input <- as.numeric(args[4])
  message("Using arguments from Command Line.")
} else {
  rep_id <- REP_ID_DEFAULT
  n_sample_input <- N_SAMPLE_DEFAULT
  l1_input <- L1_DEFAULT
  l2_input <- L2_DEFAULT
  message(sprintf("Using Defaults: Rep=%d, N=%d, L1=%.2f, L2=%.2f", rep_id, n_sample_input, l1_input, l2_input))
  message("To use command line args: Rscript script.R <ID> <N> <LAM1> <LAM2>")
}

# Validate inputs
if (is.na(rep_id) || is.na(n_sample_input) || is.na(l1_input) || is.na(l2_input)) {
  stop("Error parsing arguments. Ensure ID/N are integers and Lambdas are numeric.")
}

set.seed(rep_id)

message(sprintf("=== Starting Beta Est | Rep: %d | N: %d | L1: %.4f | L2: %.4f ===",
                rep_id, n_sample_input, l1_input, l2_input))

# ------------------------------------------------------------------------------
# 1. DATA GENERATION
# ------------------------------------------------------------------------------
generate_complex_data <- function(n_sample = 200) {
  eval_points <- seq(0, 1, length.out=100)
  basis <- create.bspline.basis(c(0,1), 20)
  nbasis <- basis$nbasis

  # Truth
  beta1_vec <- 2 * eval_points * sin(3 * pi * eval_points)
  beta1_fd  <- Data2fd(eval_points, beta1_vec, basis)

  # "Easier" Single Bump
  beta2_vec <- 2 * exp(-10 * (eval_points - 0.5)^2)
  beta2_fd  <- Data2fd(eval_points, beta2_vec, basis)

  beta_z_true <- c(1.5, -1.0)

  # Predictors
  coefs1 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- 0.6 * coefs2 + 0.4 * coefs1 # Induce correlation

  fd1 <- fd(coefs1, basis)
  fd2 <- fd(coefs2, basis)

  z_latent <- t(coefs1[1:2, ])
  Z <- z_latent + matrix(rnorm(n_sample*2, sd=0.5), n_sample, 2)

  # Response
  y_func1 <- fda::inprod(fd1, beta1_fd)
  y_func2 <- fda::inprod(fd2, beta2_fd)
  y_scalar <- Z %*% beta_z_true
  y_clean <- as.vector(y_func1 + y_func2 + y_scalar)
  y <- y_clean + rnorm(n_sample, sd=sd(y_clean)*0.05)

  # Hybrid Object
  W <- predictor_hybrid(Z, list(fd1, fd2), eval_points)

  return(list(
    W = W, y = y,
    truth = list(beta_fd = list(beta1_fd, beta2_fd), beta_z = beta_z_true)
  ))
}

# ------------------------------------------------------------------------------
# 2. RUN SIMULATION
# ------------------------------------------------------------------------------
n_comp_fixed <- 10

# Generate Data
sim <- generate_complex_data(n_sample_input)

# Pre-compute norms of true coefficients
norm_t1 <- sqrt(fda::inprod(sim$truth$beta_fd[[1]], sim$truth$beta_fd[[1]]))
norm_t2 <- sqrt(fda::inprod(sim$truth$beta_fd[[2]], sim$truth$beta_fd[[2]]))
norm_ts <- sqrt(sum(sim$truth$beta_z^2))

# Fit Model with the specific input lambdas
fit <- fit.hybridPLS(sim$W, sim$y, n_iter = n_comp_fixed, lambda = c(l1_input, l2_input))
beta_est <- fit$beta[[n_comp_fixed]]

# Calculate Errors
d1 <- minus.fd(sim$truth$beta_fd[[1]], beta_est$functional_list[[1]])
d2 <- minus.fd(sim$truth$beta_fd[[2]], beta_est$functional_list[[2]])
ds <- sim$truth$beta_z - as.vector(beta_est$Z)

# Store results
rep_results <- data.frame(
  Rep = rep_id,
  N_Sample = n_sample_input,
  Lam1 = l1_input,
  Lam2 = l2_input,
  Err_Beta1 = sqrt(fda::inprod(d1, d1)) / norm_t1,
  Err_Beta2 = sqrt(fda::inprod(d2, d2)) / norm_t2,
  Err_Scalar = sqrt(sum(ds^2)) / norm_ts
)

# ------------------------------------------------------------------------------
# 3. SAVE RESULTS (RDS)
# ------------------------------------------------------------------------------
# Filename includes Lambdas to prevent overwriting when running parameter sweeps
# Using format to avoid long floating point strings in filename
output_file <- file.path(save_dir, sprintf("beta_est_n%d_L1_%g_L2_%g_rep_%d.rds",
                                           n_sample_input, l1_input, l2_input, rep_id))
saveRDS(rep_results, output_file)

message(sprintf("Success. Saved to %s", output_file))
