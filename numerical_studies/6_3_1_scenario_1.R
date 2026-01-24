# ==============================================================================
# 6_3_1_scenario_1.R
# Section 6.3
# Scenario 1: Orthogonal Nuisance Variance (Section 6.3)
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
# 0. CONFIGURATION
# ------------------------------------------------------------------------------
# Default Configuration
REPO_DIR <- "."
OUTPUT_DIR_REL <- "./results/dgp_orthogonal"
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
set.seed(rep_id * 555)

# Increase Sample Size slightly to help separate signal from noise
n_train <- 200
n_test  <- 200
n_comp_max <- 5

message(sprintf("=== [Rep %d] Job Started. Robust Orthogonal Setting ===", rep_id))

# ------------------------------------------------------------------------------
# 1. DATA GENERATION: STRICT ORTHOGONALITY
# ------------------------------------------------------------------------------
generate_strict_orthogonal_data <- function(n_sample) {
  # Basis setup
  eval_points <- seq(0, 1, length.out=100)
  basis <- create.bspline.basis(c(0,1), 20)

  # --- STEP 1: SIGNAL & RESPONSE ---
  u_signal <- rnorm(n_sample, sd = 1)
  y_clean <- 2 * u_signal
  y <- y_clean + rnorm(n_sample, sd = sd(y_clean) * 0.05)

  # --- STEP 2: NOISE FACTOR (FORCED ORTHOGONAL TO Y) ---
  u_noise_raw <- rnorm(n_sample, sd = 1)
  mod <- lm(u_noise_raw ~ y)
  u_noise_ortho <- residuals(mod)
  u_noise <- scale(u_noise_ortho) * 5

  # --- STEP 3: FUNCTIONAL SHAPES ---
  beta1_vec <- sin(2 * pi * eval_points)
  beta2_vec <- cos(2 * pi * eval_points)
  noise1_vec <- sin(4 * pi * eval_points)
  noise2_vec <- cos(4 * pi * eval_points)

  # --- STEP 4: CONSTRUCT PREDICTORS ---
  # X1 = (Noise * Noise_Shape) + (Signal * Signal_Shape)
  X1_mat <- outer(as.vector(u_noise), noise1_vec) +
            outer(as.vector(u_signal), beta1_vec)
  X1_mat <- X1_mat + matrix(rnorm(n_sample*100, sd=0.1), n_sample, 100)
  fd1 <- Data2fd(eval_points, t(X1_mat), basis)

  X2_mat <- outer(as.vector(u_noise), noise2_vec) +
            outer(as.vector(u_signal), beta2_vec)
  X2_mat <- X2_mat + matrix(rnorm(n_sample*100, sd=0.1), n_sample, 100)
  fd2 <- Data2fd(eval_points, t(X2_mat), basis)

  # Z: Scalar Predictors
  z1 <- u_noise + rnorm(n_sample, sd=0.1)
  z2 <- u_noise + rnorm(n_sample, sd=0.1)
  z3 <- u_noise + rnorm(n_sample, sd=0.1)
  z4 <- u_noise + rnorm(n_sample, sd=0.1)
  z5 <- u_signal + rnorm(n_sample, sd=0.1)
  Z <- cbind(z1, z2, z3, z4, z5)
  colnames(Z) <- paste0("Z", 1:5)

  W <- predictor_hybrid(Z, list(fd1, fd2), eval_points)

  return(list(W = W, y = y))
}

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
subset_hybrid <- function(W, idx) {
  Z_sub <- W$Z[idx, , drop=FALSE]
  fd_list_sub <- lapply(W$functional_list, function(fd_obj) fd_obj[idx])
  predictor_hybrid(Z_sub, fd_list_sub, W$eval_point)
}

fit_hybrid_pcr_iterative <- function(W_train, y_train, W_test, y_test, n_comp_max = 10) {
  validation_rmse <- numeric(n_comp_max)

  # Storage for correlations for each component L
  cor_F1_F2 <- numeric(n_comp_max)
  cor_F1_Z  <- numeric(n_comp_max)
  cor_F2_Z  <- numeric(n_comp_max)

  # 1. Scalar PCA
  k_z_max <- min(n_comp_max, ncol(W_train$Z))
  pca_z <- prcomp(W_train$Z, center = TRUE, scale. = TRUE)
  z_tr <- pca_z$x[, 1:k_z_max, drop = FALSE]
  z_te <- predict(pca_z, newdata = W_test$Z)[, 1:k_z_max, drop = FALSE]

  # 2. Functional PCA (Robust)
  n_funcs <- length(W_train$functional_list) # Expecting 2 for this scenario
  f_tr_list <- list(); f_te_list <- list()

  for(i in 1:n_funcs) {
    pca_fd <- pca.fd(W_train$functional_list[[i]], nharm = n_comp_max)
    f_tr_list[[i]] <- pca_fd$scores

    mean_fd <- pca_fd$meanfd; harmonics <- pca_fd$harmonics
    coefs_centered <- sweep(W_test$functional_list[[i]]$coefs, 1, as.vector(mean_fd$coefs), "-")
    centered_test  <- fd(coefs_centered, W_test$functional_list[[i]]$basis)
    f_te_list[[i]] <- inprod(centered_test, harmonics)
  }

  # 3. Regression & Correlation Tracking
  for (l in 1:n_comp_max) {
    # -- A. Calculate Correlations of the l-th components --
    # Scores for functional predictor 1 and 2 at component l
    score_f1 <- f_tr_list[[1]][, l]
    score_f2 <- f_tr_list[[2]][, l]

    # Score for scalar Z at component l (handle if l > k_z_max)
    if (l <= k_z_max) {
      score_z <- z_tr[, l]
      cor_F1_Z[l] <- cor(score_f1, score_z)
      cor_F2_Z[l] <- cor(score_f2, score_z)
    } else {
      cor_F1_Z[l] <- NA
      cor_F2_Z[l] <- NA
    }

    cor_F1_F2[l] <- cor(score_f1, score_f2)

    # -- B. Regression --
    k_curr <- min(l, k_z_max)
    z_tr_c <- z_tr[, 1:k_curr, drop=FALSE]
    z_te_c <- z_te[, 1:k_curr, drop=FALSE]

    f_tr_c <- lapply(f_tr_list, function(m) m[, 1:l, drop=FALSE])
    f_te_c <- lapply(f_te_list, function(m) m[, 1:l, drop=FALSE])

    df_tr <- data.frame(z_tr_c, do.call(cbind, f_tr_c), Y = y_train)
    df_te <- data.frame(z_te_c, do.call(cbind, f_te_c))

    suppressWarnings({
      mod <- lm(Y ~ ., data = df_tr)
      pred <- predict(mod, newdata = df_te)
    })
    validation_rmse[l] <- sqrt(mean((y_test - pred)^2))
  }

  return(list(
    validation_rmse = validation_rmse,
    cor_F1_F2 = cor_F1_F2,
    cor_F1_Z = cor_F1_Z,
    cor_F2_Z = cor_F2_Z
  ))
}

tune_hybrid_pls <- function(W_train, y_train, n_folds=5, n_comp=10) {
  lambda_grid <- list(c(0,0), c(0.01,0.01), c(0.1,0.1))
  n <- length(y_train)
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  grid_results <- numeric(length(lambda_grid))

  for (g in seq_along(lambda_grid)) {
    lam <- lambda_grid[[g]]
    mses <- numeric(n_folds)
    for (k in 1:n_folds) {
      idx_val <- which(fold_ids == k); idx_tr <- which(fold_ids != k)
      w_tr <- subset_hybrid(W_train, idx_tr); y_tr <- y_train[idx_tr]
      w_val <- subset_hybrid(W_train, idx_val); y_val <- y_train[idx_val]

      tryCatch({
        fit <- fit.hybridPLS(w_tr, y_tr, n_iter=n_comp, lambda=lam,
                             validation_data = list(W_test=w_val, y_test=y_val))
        mses[k] <- min(fit$validation_rmse^2)
      }, error = function(e) mses[k] <- Inf)
    }
    grid_results[g] <- mean(mses)
  }
  return(lambda_grid[[which.min(grid_results)]])
}

# ------------------------------------------------------------------------------
# 3. MAIN LOOP
# ------------------------------------------------------------------------------

# A. Generate Data
message("   Generating Data...")
sim_data <- generate_strict_orthogonal_data(n_train + n_test)

# B. Split & Normalize
processed <- split_and_normalize.all(sim_data$W, sim_data$y, train_ratio = n_train/(n_train+n_test))
W_train <- processed$predictor_train; W_test <- processed$predictor_test
y_train <- processed$response_train;  y_test <- processed$response_test

scale_y <- sd(y_test)
sst <- sum((y_test - mean(y_test))^2)

# C. Hybrid PLS
message("   Running HybridPLS...")
best_lam <- tune_hybrid_pls(W_train, y_train, n_folds=5, n_comp=n_comp_max)
fit_pls <- fit.hybridPLS(W_train, y_train, n_iter=n_comp_max, lambda=best_lam,
                         validation_data = list(W_test = W_test, y_test = y_test))
pls_rmse <- fit_pls$validation_rmse

df_pls <- data.frame(
  Rep = rep_id, Method = "HybridPLS",
  Best_L1 = best_lam[1], Best_L2 = best_lam[2],
  Component = 1:n_comp_max,
  MSE = pls_rmse^2, Scaled_RMSE = pls_rmse / scale_y,
  R2 = 1 - ((pls_rmse^2 * n_test) / sst),
  # Fill correlations as NA for PLS
  Corr_F1_F2 = NA,
  Corr_F1_Z = NA,
  Corr_F2_Z = NA
)

# D. Hybrid PCR
message("   Running HybridPCR...")
fit_pcr <- fit_hybrid_pcr_iterative(W_train, y_train, W_test, y_test, n_comp_max)
pcr_rmse <- fit_pcr$validation_rmse

df_pcr <- data.frame(
  Rep = rep_id, Method = "HybridPCR",
  Best_L1 = NA, Best_L2 = NA,
  Component = 1:n_comp_max,
  MSE = pcr_rmse^2, Scaled_RMSE = pcr_rmse / scale_y,
  R2 = 1 - ((pcr_rmse^2 * n_test) / sst),
  # Add the correlation metrics
  Corr_F1_F2 = fit_pcr$cor_F1_F2,
  Corr_F1_Z  = fit_pcr$cor_F1_Z,
  Corr_F2_Z  = fit_pcr$cor_F2_Z
)

final_df <- rbind(df_pls, df_pcr)

# ------------------------------------------------------------------------------
# 4. SAVE TO DISK
# ------------------------------------------------------------------------------
# output_dir is already defined in Section 0


output_filename <- file.path(output_dir, sprintf("res_ortho_rep_%d.rds", rep_id))
saveRDS(final_df, file = output_filename)

message(sprintf("=== [Rep %d] Job Complete. Saved to %s ===", rep_id, output_filename))
