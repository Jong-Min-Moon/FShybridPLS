# ==============================================================================
# 6_3_2_scenario_2.R
# Section 6.3
# Scenario 2: Function-Driven Intermodal Correlation (Section 6.3)
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
# 0. SETUP
# ------------------------------------------------------------------------------
# Default Configuration
REPO_DIR <- "."
OUTPUT_DIR_REL <- "./results/dgp_usual"
REP_ID_DEFAULT <- 1

# Create directories
if (!dir.exists(file.path(REPO_DIR, "results"))) dir.create(file.path(REPO_DIR, "results"))
output_dir <- file.path(REPO_DIR, OUTPUT_DIR_REL)
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  rep_id <- as.integer(args[1])
  message(sprintf("Using Rep ID from Command Line: %d", rep_id))
} else {
  rep_id <- REP_ID_DEFAULT
  message(sprintf("Using Default Rep ID: %d", rep_id))
}

set.seed(rep_id)

message(sprintf("=== [Rep %d] Job Started. Saving to %s ===", rep_id, output_dir))

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
subset_hybrid <- function(W, idx) {
  Z_sub <- W$Z[idx, , drop=FALSE]
  fd_list_sub <- lapply(W$functional_list, function(fd_obj) fd_obj[idx])
  predictor_hybrid(Z_sub, fd_list_sub, W$eval_points)
}

# --- PCR Fitting Function (With Correlation Tracking) ---
fit_hybrid_pcr_iterative <- function(W_train, y_train, W_test, y_test, n_comp_max = 10) {

  validation_rmse <- numeric(n_comp_max)

  # Storage for correlations
  cor_F1_F2 <- numeric(n_comp_max)
  cor_F1_Z  <- numeric(n_comp_max)
  cor_F2_Z  <- numeric(n_comp_max)

  # 1. Scalar PCA (Train)
  k_z_max <- min(n_comp_max, ncol(W_train$Z))
  pca_z <- prcomp(W_train$Z, center = TRUE, scale. = TRUE)

  all_z_scores_train <- pca_z$x[, 1:k_z_max, drop = FALSE]
  all_z_scores_test  <- predict(pca_z, newdata = W_test$Z)[, 1:k_z_max, drop = FALSE]

  # 2. Functional PCA (Train)
  n_funcs <- length(W_train$functional_list)
  all_f_scores_train_list <- list()
  all_f_scores_test_list  <- list()

  for(i in 1:n_funcs) {
    # Fit FPCA on Training Data
    pca_fd <- pca.fd(W_train$functional_list[[i]], nharm = n_comp_max)
    all_f_scores_train_list[[i]] <- pca_fd$scores

    # Project Test Data (Manual Centering)
    mean_fd <- pca_fd$meanfd
    harmonics <- pca_fd$harmonics

    coefs_test <- W_test$functional_list[[i]]$coefs
    coefs_mean <- as.vector(mean_fd$coefs)

    # Sweep subtract mean
    coefs_centered <- sweep(coefs_test, 1, coefs_mean, "-")
    centered_test  <- fd(coefs_centered, W_test$functional_list[[i]]$basis)

    all_f_scores_test_list[[i]] <- inprod(centered_test, harmonics)
  }

  # 3. Iterative Regression & Correlation Calculation
  for (l in 1:n_comp_max) {

    # --- A. Calculate Correlations of the l-th components ---
    # Get the l-th score column for Func 1 and Func 2
    score_f1 <- all_f_scores_train_list[[1]][, l]
    score_f2 <- all_f_scores_train_list[[2]][, l]

    # Get l-th score for Scalar (handle if l > max scalar comps)
    if (l <= k_z_max) {
      score_z <- all_z_scores_train[, l]
      cor_F1_Z[l] <- cor(score_f1, score_z)
      cor_F2_Z[l] <- cor(score_f2, score_z)
    } else {
      cor_F1_Z[l] <- NA
      cor_F2_Z[l] <- NA
    }

    cor_F1_F2[l] <- cor(score_f1, score_f2)

    # --- B. Regression ---
    # Scalar Scores
    k_curr <- min(l, k_z_max)
    z_train_curr <- all_z_scores_train[, 1:k_curr, drop=FALSE]
    z_test_curr  <- all_z_scores_test[, 1:k_curr, drop=FALSE]

    # Functional Scores
    f_train_curr_list <- list()
    f_test_curr_list  <- list()
    for(i in 1:n_funcs) {
      f_train_curr_list[[i]] <- all_f_scores_train_list[[i]][, 1:l, drop=FALSE]
      f_test_curr_list[[i]]  <- all_f_scores_test_list[[i]][, 1:l, drop=FALSE]
    }

    # Combine
    df_train <- data.frame(z_train_curr, do.call(cbind, f_train_curr_list))
    df_test  <- data.frame(z_test_curr,  do.call(cbind, f_test_curr_list))

    df_train$Y <- y_train
    # Use suppressWarnings for rank-deficient fits if L is high
    suppressWarnings({
      model <- lm(Y ~ ., data = df_train)
      y_pred <- predict(model, newdata = df_test)
    })

    validation_rmse[l] <- sqrt(mean((y_test - y_pred)^2))
  }

  return(list(
    validation_rmse = validation_rmse,
    cor_F1_F2 = cor_F1_F2,
    cor_F1_Z = cor_F1_Z,
    cor_F2_Z = cor_F2_Z
  ))
}

# --- PLS CV Selector ---
select_best_lambda <- function(W, y, k_fold=5, comp_max=10, rep_id_log=NA) {
  lam_vals <- c(0, 0.001, 0.1, 1.0)
  grid <- expand.grid(L1 = lam_vals, L2 = lam_vals)
  n <- length(y)
  folds <- split(sample(1:n), cut(seq(1, n), breaks=k_fold, labels=FALSE))
  grid$Mean_CV_MSE <- 0

  for(g in 1:nrow(grid)) {
    l1 <- grid$L1[g]; l2 <- grid$L2[g]
    fold_mses <- numeric(k_fold)

    for(f in 1:k_fold) {
      idx_val <- folds[[f]]
      idx_tr  <- setdiff(1:n, idx_val)
      w_tr <- subset_hybrid(W, idx_tr); y_tr <- y[idx_tr]
      w_val <- subset_hybrid(W, idx_val); y_val <- y[idx_val]

      tryCatch({
        fit <- fit.hybridPLS(w_tr, y_tr, n_iter=comp_max, lambda=c(l1, l2),
                             validation_data = list(W_test = w_val, y_test = y_val))
        fold_mses[f] <- min(fit$validation_rmse^2)
      }, error = function(e) fold_mses[f] <- Inf)
    }
    grid$Mean_CV_MSE[g] <- mean(fold_mses)
  }
  best_row <- grid[which.min(grid$Mean_CV_MSE), ]
  return(c(best_row$L1, best_row$L2))
}

# ------------------------------------------------------------------------------
# 2. DATA GENERATION (ORTHOGONAL NOISE SCENARIO)
# ------------------------------------------------------------------------------
generate_torture_data <- function(n_sample) {
  eval_points <- seq(0, 1, length.out=100)
  basis <- create.bspline.basis(c(0,1), 20)
  nbasis <- basis$nbasis

  beta1_vec <- 2 * eval_points * sin(3 * pi * eval_points)
  beta1_fd  <- Data2fd(eval_points, beta1_vec, basis)
  beta2_vec <- 2 * exp(-10 * (eval_points - 0.5)^2)
  beta2_fd  <- Data2fd(eval_points, beta2_vec, basis)
  beta_z_true <- c(1.5, -1.0, 1.5, -1.0, 1.5, -1.0)/5

  coefs1 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- 0.6 * coefs2 + 0.4 * coefs1

  fd1 <- fd(coefs1, basis)
  fd2 <- fd(coefs2, basis)

  z_latent <- t(coefs1[1:6, ])
  Z <- z_latent + matrix(rnorm(n_sample*2, sd=0.5), n_sample, 6)

  y_func1 <- inprod(fd1, beta1_fd)
  y_func2 <- inprod(fd2, beta2_fd)
  y_scalar <- Z %*% beta_z_true
  y_clean <- as.vector(y_func1 + y_func2 + y_scalar)
  y <- y_clean + rnorm(n_sample, sd=sd(y_clean)*0.05)

  W <- predictor_hybrid(Z, list(fd1, fd2), eval_points)
  return(list(W = W, y = y))
}

# ------------------------------------------------------------------------------
# 3. MAIN SIMULATION LOOP
# ------------------------------------------------------------------------------
n_train <- 100
n_test  <- 100
n_comp_max <- 5

message(sprintf("   [Rep %d] Generating Data...", rep_id))
sim_data <- generate_torture_data(n_train + n_test)
W_train <- subset_hybrid(sim_data$W, 1:n_train)
y_train <- sim_data$y[1:n_train]
W_test  <- subset_hybrid(sim_data$W, (n_train+1):(n_train+n_test))
y_test  <- sim_data$y[(n_train+1):(n_train+n_test)]

# Pre-calculate stats for R2
scale_y <- sd(y_test)
sst <- sum((y_test - mean(y_test))^2)

# --- METHOD 1: Hybrid PLS ---
message(sprintf("   [Rep %d] Running HybridPLS...", rep_id))
best_lam <- select_best_lambda(W_train, y_train, k_fold=5, comp_max=n_comp_max, rep_id_log=rep_id)
fit_pls <- fit.hybridPLS(W_train, y_train, n_iter=n_comp_max, lambda=best_lam,
                         validation_data = list(W_test = W_test, y_test = y_test))

pls_rmse <- fit_pls$validation_rmse
pls_mse  <- pls_rmse^2
pls_r2   <- 1 - ((pls_mse * n_test) / sst)

df_pls <- data.frame(
  Rep = rep_id,
  Method = "HybridPLS",
  Best_L1 = best_lam[1],
  Best_L2 = best_lam[2],
  Component = 1:n_comp_max,
  MSE = pls_mse,
  Scaled_RMSE = pls_rmse / scale_y,
  R2 = pls_r2,
  # NA for Correlations
  Corr_F1_F2 = NA,
  Corr_F1_Z = NA,
  Corr_F2_Z = NA
)

# --- METHOD 2: Hybrid PCR ---
message(sprintf("   [Rep %d] Running HybridPCR...", rep_id))
fit_pcr <- fit_hybrid_pcr_iterative(W_train, y_train, W_test, y_test, n_comp_max=n_comp_max)

pcr_rmse <- fit_pcr$validation_rmse
pcr_mse  <- pcr_rmse^2
pcr_r2   <- 1 - ((pcr_mse * n_test) / sst)

df_pcr <- data.frame(
  Rep = rep_id,
  Method = "HybridPCR",
  Best_L1 = NA,
  Best_L2 = NA,
  Component = 1:n_comp_max,
  MSE = pcr_mse,
  Scaled_RMSE = pcr_rmse / scale_y,
  R2 = pcr_r2,
  # Save Correlations
  Corr_F1_F2 = fit_pcr$cor_F1_F2,
  Corr_F1_Z  = fit_pcr$cor_F1_Z,
  Corr_F2_Z  = fit_pcr$cor_F2_Z
)

# --- SAVE ALL ---
results_log <- rbind(df_pls, df_pcr)

# Create filename
output_filename <- file.path(output_dir, sprintf("res_ortho_rep_%d.rds", rep_id))

# Save
saveRDS(results_log, file = output_filename)
message(sprintf("=== [Rep %d] Success: Saved to %s ===", rep_id, output_filename))
