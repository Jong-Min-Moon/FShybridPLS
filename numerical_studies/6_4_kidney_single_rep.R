# ==============================================================================
# 6_4_kidney_single_rep.R
# Section 6.4
# Run ONE replication of HybridPLS vs HybridPCR vs OLS vs PFR
# kidney data file is not provided
# ==============================================================================
# Installation of FSHybridPLS package
# This package is generated using litr, a literate programming tool developed by Jacob Bien and Patrick Vossler (https://jacobbien.github.io/litr-project/). I used litr to build the package and you don't have to install litr just to use the package. To install the package:
#
# devtools::install("your_downloaded_directory/FSHybridPLS")

library(fda)
library(dplyr)
library(tidyr)
library(refund)      # For PFR
library(FSHybridPLS)


# ------------------------------------------------------------------------------
# 0. PARSE ARGUMENTS
# ------------------------------------------------------------------------------
# Default Configuration
REPO_DIR <- "."
SAVE_DIR_REL <- "./results/real_data"
DATA_DIR_REL <- "../data"  # Assuming data is in a sibling directory called 'data'

LAM1_DEFAULT <- 0.1
LAM2_DEFAULT <- 0.01
REP_ID_DEFAULT <- 1

# Create 'results' directory if it doesn't exist
if (!dir.exists(file.path(REPO_DIR, "results"))) dir.create(file.path(REPO_DIR, "results"))
object_dir <- file.path(REPO_DIR, SAVE_DIR_REL)
if(!dir.exists(object_dir)) dir.create(object_dir, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 3) {
  lam1   <- as.numeric(args[1])
  lam2   <- as.numeric(args[2])
  rep_id <- as.integer(args[3])
  message("Using arguments from Command Line.")
} else {
  lam1   <- LAM1_DEFAULT
  lam2   <- LAM2_DEFAULT
  rep_id <- REP_ID_DEFAULT
  message(sprintf("Using Defaults: L1=%.2f, L2=%.2f, Rep=%d", lam1, lam2, rep_id))
}

lambda_vec <- c(lam1, lam2)

message(sprintf("=== Job Started. Rep: %d | L1: %.4f | L2: %.4f ===", rep_id, lam1, lam2))

# ------------------------------------------------------------------------------
# 1. HELPER FUNCTIONS
# ------------------------------------------------------------------------------
subset_hybrid <- function(W, idx) {
  Z_sub <- W$Z[idx, , drop=FALSE]
  fd_list_sub <- lapply(W$functional_list, function(fd_obj) fd_obj[idx])
  predictor_hybrid(Z_sub, fd_list_sub, W$eval_point)
}

# --- PCR Wrapper ---
fit_hybrid_pcr_iterative_objects <- function(W_train, y_train, W_test, y_test, n_comp_max = 10) {
  validation_rmse <- numeric(n_comp_max)

  cor_F1_F2 <- numeric(n_comp_max)
  cor_F1_Z  <- numeric(n_comp_max)
  cor_F2_Z  <- numeric(n_comp_max)

  k_z_max <- min(n_comp_max, ncol(W_train$Z))
  pca_z <- prcomp(W_train$Z, center = TRUE, scale. = TRUE)

  z_tr <- pca_z$x[, 1:k_z_max, drop = FALSE]
  z_te <- predict(pca_z, newdata = W_test$Z)[, 1:k_z_max, drop = FALSE]

  n_funcs <- length(W_train$functional_list)
  f_tr_list <- list(); f_te_list <- list()
  pca_fd_list <- list()

  for(i in 1:n_funcs) {
    pca_fd <- pca.fd(W_train$functional_list[[i]], nharm = n_comp_max)
    pca_fd_list[[i]] <- pca_fd

    f_tr_list[[i]] <- pca_fd$scores

    mean_fd <- pca_fd$meanfd; harmonics <- pca_fd$harmonics
    coefs_centered <- sweep(W_test$functional_list[[i]]$coefs, 1, as.vector(mean_fd$coefs), "-")
    centered_test  <- fd(coefs_centered, W_test$functional_list[[i]]$basis)
    f_te_list[[i]] <- inprod(centered_test, harmonics)
  }

  for (l in 1:n_comp_max) {
    score_f1 <- f_tr_list[[1]][, l]
    score_f2 <- f_tr_list[[2]][, l]

    if (l <= k_z_max) {
      score_z <- z_tr[, l]
      cor_F1_Z[l] <- cor(score_f1, score_z)
      cor_F2_Z[l] <- cor(score_f2, score_z)
    } else {
      cor_F1_Z[l] <- NA
      cor_F2_Z[l] <- NA
    }
    cor_F1_F2[l] <- cor(score_f1, score_f2)

    k_curr <- min(l, k_z_max)
    df_tr <- data.frame(z_tr[, 1:k_curr, drop=FALSE],
                        do.call(cbind, lapply(f_tr_list, function(m) m[, 1:l, drop=FALSE])),
                        Y = y_train)
    df_te <- data.frame(z_te[, 1:k_curr, drop=FALSE],
                        do.call(cbind, lapply(f_te_list, function(m) m[, 1:l, drop=FALSE])))

    suppressWarnings({ mod <- lm(Y ~ ., data = df_tr); pred <- predict(mod, newdata = df_te) })
    validation_rmse[l] <- sqrt(mean((y_test - pred)^2))
  }

  return(list(
    validation_rmse = validation_rmse,
    pca_z_object = pca_z,
    pca_fd_objects = pca_fd_list,
    correlations = data.frame(Comp=1:n_comp_max, F1_F2=cor_F1_F2, F1_Z=cor_F1_Z, F2_Z=cor_F2_Z)
  ))
}

# --- PFR Wrapper (Refund) ---
# Now returns a list with RMSE and the Fit Object
fit_pfr_legacy <- function(W_train, y_train, W_test, y_test) {
  eval_pts <- W_train$eval_point

  Z_train_df <- as.data.frame(W_train$Z)
  colnames(Z_train_df) <- paste0("Z", 1:ncol(Z_train_df))

  Z_test_df  <- as.data.frame(W_test$Z)
  colnames(Z_test_df) <- paste0("Z", 1:ncol(Z_test_df))

  df_train <- cbind(response = y_train, Z_train_df)
  df_train$F1 <- t(eval.fd(eval_pts, W_train$functional_list[[1]]))
  df_train$F2 <- t(eval.fd(eval_pts, W_train$functional_list[[2]]))

  df_test <- cbind(response = y_test, Z_test_df)
  df_test$F1 <- t(eval.fd(eval_pts, W_test$functional_list[[1]]))
  df_test$F2 <- t(eval.fd(eval_pts, W_test$functional_list[[2]]))

  scalar_formula_part <- paste(colnames(Z_train_df), collapse = " + ")

  f_str <- paste0(
    "response ~ ", scalar_formula_part, " + ",
    "lf(F1, argvals = eval_pts, presmooth = 'fpca.sc') + ",
    "lf(F2, argvals = eval_pts, presmooth = 'fpca.sc' )"
  )

  tryCatch({
    fit <- pfr(as.formula(f_str), data = df_train, method="GCV.Cp", gamma=1.2)
    pred <- predict(fit, newdata = df_test)
    rmse <- sqrt(mean((y_test - pred)^2))
    return(list(rmse = rmse, fit_obj = fit))
  }, error = function(e) {
    message("PFR Error: ", e$message)
    return(list(rmse = NA, fit_obj = NULL))
  })
}

# ------------------------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------------------------
# Using relative data path
data_path <- file.path(REPO_DIR, DATA_DIR_REL, "renogram_data.csv")
if(!file.exists(data_path)) stop("renogram_data.csv not found!")

raw_df <- read.csv(data_path)
kidney_data <- load_and_preprocess_kidney_data(raw_df, n_basis = 20)
W_full <- kidney_data$W
y_full <- kidney_data$y

# ------------------------------------------------------------------------------
# 3. RUN SIMULATION
# ------------------------------------------------------------------------------
n_comp_max      <- 16
subsample_ratio <- 0.40
n_total <- length(y_full)
n_sub   <- floor(n_total * subsample_ratio)

set.seed(rep_id * 123)

idx_sub <- sample(1:n_total, n_sub, replace = FALSE)
W_working <- subset_hybrid(W_full, idx_sub)
y_working <- y_full[idx_sub]

processed <- split_and_normalize.all(W_working, y_working, train_ratio = 0.7)
W_train <- processed$predictor_train; W_test <- processed$predictor_test
y_train <- processed$response_train;  y_test <- processed$response_test

y_range <- max(y_test) - min(y_test); if(y_range == 0) y_range <- 1

# --- Method 1: Hybrid PLS (Iterative) ---
message("Running Hybrid PLS...")
pls_fit <- fit.hybridPLS(W_train, y_train, n_iter = n_comp_max, lambda = lambda_vec,
                         validation_data = list(W_test = W_test, y_test = y_test))
rmse_pls <- pls_fit$validation_rmse

# --- Method 2: Hybrid PCR (Iterative) ---
message("Running Hybrid PCR...")
pcr_out <- fit_hybrid_pcr_iterative_objects(W_train, y_train, W_test, y_test, n_comp_max)
rmse_pcr <- pcr_out$validation_rmse

# --- Method 3: Scalar OLS (Run Once) ---
message("Running Scalar OLS...")
Z_train_df <- as.data.frame(W_train$Z)
Z_test_df  <- as.data.frame(W_test$Z)
df_ols_train <- cbind(Z_train_df, Y = y_train)

fit_ols <- lm(Y ~ ., data = df_ols_train)
ols_coefs <- coef(fit_ols)
pred_ols <- predict(fit_ols, newdata = Z_test_df)
rmse_ols <- sqrt(mean((y_test - pred_ols)^2))

# --- Method 4: PFR (Run Once) ---
message("Running PFR...")
pfr_res <- fit_pfr_legacy(W_train, y_train, W_test, y_test)
rmse_pfr <- pfr_res$rmse

# ------------------------------------------------------------------------------
# 4. COMPILE RESULTS
# ------------------------------------------------------------------------------
# Hybrid PLS Results
df_pls <- data.frame(Rep = rep_id, Method = "HybridPLS", L1 = lam1, L2 = lam2,
                     Component = 1:n_comp_max, NRMSE = rmse_pls / y_range,
                     Corr_F1_F2=NA, Corr_F1_Z=NA, Corr_F2_Z=NA)

# Hybrid PCR Results
df_pcr <- data.frame(Rep = rep_id, Method = "HybridPCR", L1 = lam1, L2 = lam2,
                     Component = 1:n_comp_max, NRMSE = rmse_pcr / y_range,
                     Corr_F1_F2=pcr_out$correlations$F1_F2,
                     Corr_F1_Z =pcr_out$correlations$F1_Z,
                     Corr_F2_Z =pcr_out$correlations$F2_Z)

# Scalar OLS Results (Replicated for plotting convenience)
df_ols <- data.frame(Rep = rep_id, Method = "ScalarOLS", L1 = lam1, L2 = lam2,
                     Component = 1:n_comp_max, NRMSE = rmse_ols / y_range,
                     Corr_F1_F2=NA, Corr_F1_Z=NA, Corr_F2_Z=NA)

# PFR Results (Replicated for plotting convenience)
df_pfr <- data.frame(Rep = rep_id, Method = "PFR", L1 = lam1, L2 = lam2,
                     Component = 1:n_comp_max, NRMSE = rmse_pfr / y_range,
                     Corr_F1_F2=NA, Corr_F1_Z=NA, Corr_F2_Z=NA)

numeric_results_df <- rbind(df_pls, df_pcr, df_ols, df_pfr)

# ------------------------------------------------------------------------------
# 5. SAVE ALL TO RDS
# ------------------------------------------------------------------------------
saved_data <- list(
  # Metadata
  rep_id = rep_id,
  lambda = lambda_vec,

  # A. The Main Numeric Results (Accuracy & Correlations)
  numeric_results = numeric_results_df,

  # B. Model Objects (for later detailed analysis)
  pls_fit_object   = pls_fit,
  pca_z_object     = pcr_out$pca_z_object,
  pca_fd_objects   = pcr_out$pca_fd_objects,
  ols_coefficients = ols_coefs,
  pfr_fit_object   = pfr_res$fit_obj  # Save PFR object
)

# Filename: results_rep_X.rds
output_filename <- file.path(object_dir, sprintf("results_rep_%d.rds", rep_id))
saveRDS(saved_data, file = output_filename)

message(sprintf("=== Replication %d Complete. All data saved to %s ===", rep_id, output_filename))
