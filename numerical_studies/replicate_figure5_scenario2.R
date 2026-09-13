# Figure 5 (fig:simul_2 / synthetic_usual.pdf) — Scenario 2
# Old one-rep driver: numerical_studies/6_3_2_scenario_2.R
# Old aggregated results:
#   development/.../6_3_simul_pred/synthetic_usual/result.csv
# Plot: synthetic_usual/plot.R (and/or plot/)
#
# Paper: n=200 (100/100), 100 MC, Hybrid PLS vs PCR scaled RMSE + PCR |cor|
# Note: script uses sequential split (not split_and_normalize) and 4x4 lambda CV.
suppressPackageStartupMessages({
  library(FSHybridPLS)
  library(fda)
})

n_mc <- 50L
n_train <- 100L
n_test <- 100L
n_comp_max <- 5L
out_dir <- "C:/Users/jongmin/Downloads/FShybridPLS-main/FShybridPLS-main/numerical_studies"
old_csv <- "D:/GitHub/jongmin/FSHybridPLS_actual_code/development/simul_new/6_3_simul_pred/synthetic_usual/result.csv"

cat("=== Figure 5 / Scenario 2 (function-driven cross-modality) ===\n")
cat("Package:", as.character(packageVersion("FSHybridPLS")), "\n")
cat("MC =", n_mc, "| n_train =", n_train, "| n_test =", n_test, "| L_max =", n_comp_max, "\n\n")

subset_hybrid <- function(W, idx) {
  predictor_hybrid(
    W$Z[idx, , drop = FALSE],
    lapply(W$functional_list, function(fd_obj) fd_obj[idx]),
    W$eval_point
  )
}

fit_hybrid_pcr_iterative <- function(W_train, y_train, W_test, y_test, n_comp_max = 10) {
  validation_rmse <- numeric(n_comp_max)
  cor_F1_F2 <- numeric(n_comp_max)
  cor_F1_Z <- numeric(n_comp_max)
  cor_F2_Z <- numeric(n_comp_max)
  k_z_max <- min(n_comp_max, ncol(W_train$Z))
  pca_z <- prcomp(W_train$Z, center = TRUE, scale. = TRUE)
  all_z_scores_train <- pca_z$x[, 1:k_z_max, drop = FALSE]
  all_z_scores_test <- predict(pca_z, newdata = W_test$Z)[, 1:k_z_max, drop = FALSE]
  n_funcs <- length(W_train$functional_list)
  all_f_scores_train_list <- list()
  all_f_scores_test_list <- list()
  for (i in seq_len(n_funcs)) {
    pca_fd <- pca.fd(W_train$functional_list[[i]], nharm = n_comp_max)
    all_f_scores_train_list[[i]] <- pca_fd$scores
    mean_fd <- pca_fd$meanfd
    harmonics <- pca_fd$harmonics
    coefs_centered <- sweep(
      W_test$functional_list[[i]]$coefs, 1, as.vector(mean_fd$coefs), "-"
    )
    centered_test <- fd(coefs_centered, W_test$functional_list[[i]]$basis)
    all_f_scores_test_list[[i]] <- inprod(centered_test, harmonics)
  }
  for (l in seq_len(n_comp_max)) {
    score_f1 <- all_f_scores_train_list[[1]][, l]
    score_f2 <- all_f_scores_train_list[[2]][, l]
    if (l <= k_z_max) {
      score_z <- all_z_scores_train[, l]
      cor_F1_Z[l] <- cor(score_f1, score_z)
      cor_F2_Z[l] <- cor(score_f2, score_z)
    } else {
      cor_F1_Z[l] <- NA
      cor_F2_Z[l] <- NA
    }
    cor_F1_F2[l] <- cor(score_f1, score_f2)
    k_curr <- min(l, k_z_max)
    df_train <- data.frame(
      all_z_scores_train[, 1:k_curr, drop = FALSE],
      do.call(cbind, lapply(all_f_scores_train_list, function(m) m[, 1:l, drop = FALSE])),
      Y = y_train
    )
    df_test <- data.frame(
      all_z_scores_test[, 1:k_curr, drop = FALSE],
      do.call(cbind, lapply(all_f_scores_test_list, function(m) m[, 1:l, drop = FALSE]))
    )
    suppressWarnings({
      model <- lm(Y ~ ., data = df_train)
      y_pred <- predict(model, newdata = df_test)
    })
    validation_rmse[l] <- sqrt(mean((y_test - y_pred)^2))
  }
  list(
    validation_rmse = validation_rmse,
    cor_F1_F2 = cor_F1_F2,
    cor_F1_Z = cor_F1_Z,
    cor_F2_Z = cor_F2_Z
  )
}

select_best_lambda <- function(W, y, k_fold = 5, comp_max = 10) {
  lam_vals <- c(0, 0.001, 0.1, 1.0)
  grid <- expand.grid(L1 = lam_vals, L2 = lam_vals)
  n <- length(y)
  folds <- split(sample(1:n), cut(seq(1, n), breaks = k_fold, labels = FALSE))
  grid$Mean_CV_MSE <- 0
  for (g in seq_len(nrow(grid))) {
    l1 <- grid$L1[g]
    l2 <- grid$L2[g]
    fold_mses <- numeric(k_fold)
    for (f in seq_len(k_fold)) {
      idx_val <- folds[[f]]
      idx_tr <- setdiff(seq_len(n), idx_val)
      w_tr <- subset_hybrid(W, idx_tr)
      y_tr <- y[idx_tr]
      w_val <- subset_hybrid(W, idx_val)
      y_val <- y[idx_val]
      tryCatch({
        fit <- fit.hybridPLS(
          w_tr, y_tr, n_iter = comp_max, lambda = c(l1, l2),
          validation_data = list(W_test = w_val, y_test = y_val)
        )
        fold_mses[f] <- min(fit$validation_rmse^2)
      }, error = function(e) fold_mses[f] <<- Inf)
    }
    grid$Mean_CV_MSE[g] <- mean(fold_mses)
  }
  best_row <- grid[which.min(grid$Mean_CV_MSE), ]
  c(best_row$L1, best_row$L2)
}

# Matches 6_3_2_scenario_2.R (including Z noise size quirk: n_sample*2 recycled into nx6)
generate_torture_data <- function(n_sample) {
  eval_points <- seq(0, 1, length.out = 100)
  basis <- create.bspline.basis(c(0, 1), 20)
  nbasis <- basis$nbasis
  beta1_fd <- Data2fd(eval_points, 2 * eval_points * sin(3 * pi * eval_points), basis)
  beta2_fd <- Data2fd(eval_points, 2 * exp(-10 * (eval_points - 0.5)^2), basis)
  beta_z_true <- c(1.5, -1.0, 1.5, -1.0, 1.5, -1.0) / 5
  coefs1 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- 0.6 * coefs2 + 0.4 * coefs1
  fd1 <- fd(coefs1, basis)
  fd2 <- fd(coefs2, basis)
  z_latent <- t(coefs1[1:6, ])
  Z <- z_latent + matrix(rnorm(n_sample * 2, sd = 0.5), n_sample, 6)
  y_clean <- as.vector(inprod(fd1, beta1_fd) + inprod(fd2, beta2_fd) + Z %*% beta_z_true)
  y <- y_clean + rnorm(n_sample, sd = sd(y_clean) * 0.05)
  list(W = predictor_hybrid(Z, list(fd1, fd2), eval_points), y = y)
}

one_rep <- function(rep_id) {
  set.seed(rep_id)
  sim_data <- generate_torture_data(n_train + n_test)
  W_train <- subset_hybrid(sim_data$W, seq_len(n_train))
  y_train <- sim_data$y[seq_len(n_train)]
  W_test <- subset_hybrid(sim_data$W, (n_train + 1):(n_train + n_test))
  y_test <- sim_data$y[(n_train + 1):(n_train + n_test)]
  scale_y <- sd(y_test)
  sst <- sum((y_test - mean(y_test))^2)

  best_lam <- select_best_lambda(W_train, y_train, k_fold = 5, comp_max = n_comp_max)
  fit_pls <- fit.hybridPLS(
    W_train, y_train, n_iter = n_comp_max, lambda = best_lam,
    validation_data = list(W_test = W_test, y_test = y_test)
  )
  pls_rmse <- fit_pls$validation_rmse
  df_pls <- data.frame(
    Rep = rep_id, Method = "HybridPLS",
    Best_L1 = best_lam[1], Best_L2 = best_lam[2],
    Component = seq_len(n_comp_max),
    MSE = pls_rmse^2, Scaled_RMSE = pls_rmse / scale_y,
    R2 = 1 - ((pls_rmse^2 * n_test) / sst),
    Corr_F1_F2 = NA_real_, Corr_F1_Z = NA_real_, Corr_F2_Z = NA_real_
  )

  fit_pcr <- fit_hybrid_pcr_iterative(
    W_train, y_train, W_test, y_test, n_comp_max = n_comp_max
  )
  pcr_rmse <- fit_pcr$validation_rmse
  df_pcr <- data.frame(
    Rep = rep_id, Method = "HybridPCR",
    Best_L1 = NA_real_, Best_L2 = NA_real_,
    Component = seq_len(n_comp_max),
    MSE = pcr_rmse^2, Scaled_RMSE = pcr_rmse / scale_y,
    R2 = 1 - ((pcr_rmse^2 * n_test) / sst),
    Corr_F1_F2 = fit_pcr$cor_F1_F2,
    Corr_F1_Z = fit_pcr$cor_F1_Z,
    Corr_F2_Z = fit_pcr$cor_F2_Z
  )
  rbind(df_pls, df_pcr)
}

rows <- vector("list", n_mc)
for (r in seq_len(n_mc)) {
  rows[[r]] <- one_rep(r)
  if (r %% 5L == 0L || r == 1L) {
    cat(sprintf("  finished rep %d / %d\n", r, n_mc))
  }
}
ours <- do.call(rbind, rows)
write.csv(ours, file.path(out_dir, "replicate_figure5_scenario2_results.csv"), row.names = FALSE)
saveRDS(ours, file.path(out_dir, "replicate_figure5_scenario2_results.rds"))

summarize_rmse <- function(df) {
  do.call(rbind, lapply(split(df, list(df$Method, df$Component), drop = TRUE), function(s) {
    x <- as.numeric(s$Scaled_RMSE)
    x <- x[is.finite(x)]
    data.frame(
      Method = s$Method[1],
      Component = s$Component[1],
      n = length(x),
      mean = mean(x),
      median = median(x),
      se = sd(x) / sqrt(length(x)),
      stringsAsFactors = FALSE
    )
  }))
}

summarize_corr <- function(df) {
  pcr <- df[df$Method == "HybridPCR", ]
  do.call(rbind, lapply(split(pcr, pcr$Component), function(s) {
    data.frame(
      Component = s$Component[1],
      mean_abs_F1F2 = mean(abs(as.numeric(s$Corr_F1_F2)), na.rm = TRUE),
      mean_abs_F1Z = mean(abs(as.numeric(s$Corr_F1_Z)), na.rm = TRUE),
      mean_abs_F2Z = mean(abs(as.numeric(s$Corr_F2_Z)), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

old <- read.csv(old_csv, stringsAsFactors = FALSE)
old$Scaled_RMSE <- as.numeric(old$Scaled_RMSE)
old$Corr_F1_F2 <- as.numeric(old$Corr_F1_F2)
old$Corr_F1_Z <- as.numeric(old$Corr_F1_Z)
old$Corr_F2_Z <- as.numeric(old$Corr_F2_Z)

cat("\n=== Scaled RMSE means (old vs ours) ===\n")
sum_old <- summarize_rmse(old)
sum_new <- summarize_rmse(ours)
cmp <- merge(sum_old, sum_new, by = c("Method", "Component"), suffixes = c("_old", "_new"))
cmp$d_mean <- cmp$mean_new - cmp$mean_old
cmp$rel <- abs(cmp$d_mean) / pmax(cmp$mean_old, 1e-12)
print(cmp[order(cmp$Component, cmp$Method),
          c("Method", "Component", "n_old", "n_new", "mean_old", "mean_new",
            "d_mean", "rel", "median_old", "median_new")],
      digits = 4, row.names = FALSE)

cat("\n=== PCR |cor| means (panel b) ===\n")
c_old <- summarize_corr(old)
c_new <- summarize_corr(ours)
cc <- merge(c_old, c_new, by = "Component", suffixes = c("_old", "_new"))
print(cc[order(cc$Component), ], digits = 4, row.names = FALSE)

pls1_old <- sum_old$mean[sum_old$Method == "HybridPLS" & sum_old$Component == 1]
pcr1_old <- sum_old$mean[sum_old$Method == "HybridPCR" & sum_old$Component == 1]
pls1_new <- sum_new$mean[sum_new$Method == "HybridPLS" & sum_new$Component == 1]
pcr1_new <- sum_new$mean[sum_new$Method == "HybridPCR" & sum_new$Component == 1]

cat("\n=== Checks ===\n")
cat(sprintf("Paper text ~ PLS 0.27 / PCR 0.74 at L=1\n"))
cat(sprintf("Old CSV:  PLS %.3f / PCR %.3f\n", pls1_old, pcr1_old))
cat(sprintf("Ours:     PLS %.3f / PCR %.3f\n", pls1_new, pcr1_new))
cat(sprintf("PLS << PCR at L=1: %s\n", pls1_new < 0.6 * pcr1_new))

pls2 <- sum_new$mean[sum_new$Method == "HybridPLS" & sum_new$Component == 2]
pcr2 <- sum_new$mean[sum_new$Method == "HybridPCR" & sum_new$Component == 2]
cat(sprintf("Both improve by L=2: PLS=%.3f PCR=%.3f\n", pls2, pcr2))

c1 <- c_new[c_new$Component == 1, ]
# Scenario 2: strong F-Z correlations expected
high_fz <- (c1$mean_abs_F1Z > 0.5) || (c1$mean_abs_F2Z > 0.5)
cat(sprintf("PCR cross-modality |cor| elevated at L=1: %s (F1Z=%.3f F2Z=%.3f F1F2=%.3f)\n",
            high_fz, c1$mean_abs_F1Z, c1$mean_abs_F2Z, c1$mean_abs_F1F2))

match_l1 <- abs(pls1_new - pls1_old) / pls1_old < 0.25 &&
  abs(pcr1_new - pcr1_old) / pcr1_old < 0.25
overall <- match_l1 && (pls1_new < 0.6 * pcr1_new) && high_fz
cat(sprintf("\nFIGURE 5 OVERALL vs old result.csv: %s\n",
            if (overall) "PASS" else "REVIEW"))
cat("Saved replicate_figure5_scenario2_results.{rds,csv}\n")
