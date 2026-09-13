# Figure 4 (fig:simul_1 / synthetic_latent.pdf) — Scenario 1
# Old one-rep driver: numerical_studies/6_3_1_scenario_1.R
# Old aggregated results + plot data:
#   development/.../6_3_simul_pred/synthetic_orthogonal/result.csv
#   plot from synthetic_orthogonal/plot/plot.R
#
# Paper: n=400 (200/200 train/test), 100 MC, Hybrid PLS vs PCR scaled RMSE + PCR |cor|
suppressPackageStartupMessages({
  library(FSHybridPLS)
  library(fda)
})

n_mc <- 50L
n_train <- 200L
n_test <- 200L
n_comp_max <- 5L
out_dir <- "C:/Users/jongmin/Downloads/FShybridPLS-main/FShybridPLS-main/numerical_studies"
old_csv <- "D:/GitHub/jongmin/FSHybridPLS_actual_code/development/simul_new/6_3_simul_pred/synthetic_orthogonal/result.csv"

cat("=== Figure 4 / Scenario 1 (orthogonal nuisance) ===\n")
cat("Package:", as.character(packageVersion("FSHybridPLS")), "\n")
cat("MC =", n_mc, "| n_train =", n_train, "| n_test =", n_test, "| L_max =", n_comp_max, "\n\n")

generate_strict_orthogonal_data <- function(n_sample) {
  eval_points <- seq(0, 1, length.out = 100)
  basis <- create.bspline.basis(c(0, 1), 20)
  u_signal <- rnorm(n_sample, sd = 1)
  y_clean <- 2 * u_signal
  y <- y_clean + rnorm(n_sample, sd = sd(y_clean) * 0.05)
  u_noise_raw <- rnorm(n_sample, sd = 1)
  u_noise <- scale(residuals(lm(u_noise_raw ~ y))) * 5
  beta1_vec <- sin(2 * pi * eval_points)
  beta2_vec <- cos(2 * pi * eval_points)
  noise1_vec <- sin(4 * pi * eval_points)
  noise2_vec <- cos(4 * pi * eval_points)
  X1_mat <- outer(as.vector(u_noise), noise1_vec) +
    outer(as.vector(u_signal), beta1_vec) +
    matrix(rnorm(n_sample * 100, sd = 0.1), n_sample, 100)
  X2_mat <- outer(as.vector(u_noise), noise2_vec) +
    outer(as.vector(u_signal), beta2_vec) +
    matrix(rnorm(n_sample * 100, sd = 0.1), n_sample, 100)
  fd1 <- Data2fd(eval_points, t(X1_mat), basis)
  fd2 <- Data2fd(eval_points, t(X2_mat), basis)
  z1 <- u_noise + rnorm(n_sample, sd = 0.1)
  z2 <- u_noise + rnorm(n_sample, sd = 0.1)
  z3 <- u_noise + rnorm(n_sample, sd = 0.1)
  z4 <- u_noise + rnorm(n_sample, sd = 0.1)
  z5 <- u_signal + rnorm(n_sample, sd = 0.1)
  Z <- cbind(z1, z2, z3, z4, z5)
  colnames(Z) <- paste0("Z", 1:5)
  list(W = predictor_hybrid(Z, list(fd1, fd2), eval_points), y = y)
}

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
  z_tr <- pca_z$x[, 1:k_z_max, drop = FALSE]
  z_te <- predict(pca_z, newdata = W_test$Z)[, 1:k_z_max, drop = FALSE]
  n_funcs <- length(W_train$functional_list)
  f_tr_list <- list()
  f_te_list <- list()
  for (i in seq_len(n_funcs)) {
    pca_fd <- pca.fd(W_train$functional_list[[i]], nharm = n_comp_max)
    f_tr_list[[i]] <- pca_fd$scores
    mean_fd <- pca_fd$meanfd
    harmonics <- pca_fd$harmonics
    coefs_centered <- sweep(
      W_test$functional_list[[i]]$coefs, 1, as.vector(mean_fd$coefs), "-"
    )
    centered_test <- fd(coefs_centered, W_test$functional_list[[i]]$basis)
    f_te_list[[i]] <- inprod(centered_test, harmonics)
  }
  for (l in seq_len(n_comp_max)) {
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
    df_tr <- data.frame(
      z_tr[, 1:k_curr, drop = FALSE],
      do.call(cbind, lapply(f_tr_list, function(m) m[, 1:l, drop = FALSE])),
      Y = y_train
    )
    df_te <- data.frame(
      z_te[, 1:k_curr, drop = FALSE],
      do.call(cbind, lapply(f_te_list, function(m) m[, 1:l, drop = FALSE]))
    )
    suppressWarnings({
      mod <- lm(Y ~ ., data = df_tr)
      pred <- predict(mod, newdata = df_te)
    })
    validation_rmse[l] <- sqrt(mean((y_test - pred)^2))
  }
  list(
    validation_rmse = validation_rmse,
    cor_F1_F2 = cor_F1_F2,
    cor_F1_Z = cor_F1_Z,
    cor_F2_Z = cor_F2_Z
  )
}

tune_hybrid_pls <- function(W_train, y_train, n_folds = 5, n_comp = 10) {
  lambda_grid <- list(c(0, 0), c(0.01, 0.01), c(0.1, 0.1))
  n <- length(y_train)
  fold_ids <- sample(rep(1:n_folds, length.out = n))
  grid_results <- numeric(length(lambda_grid))
  for (g in seq_along(lambda_grid)) {
    lam <- lambda_grid[[g]]
    mses <- numeric(n_folds)
    for (k in seq_len(n_folds)) {
      idx_val <- which(fold_ids == k)
      idx_tr <- which(fold_ids != k)
      w_tr <- subset_hybrid(W_train, idx_tr)
      y_tr <- y_train[idx_tr]
      w_val <- subset_hybrid(W_train, idx_val)
      y_val <- y_train[idx_val]
      tryCatch({
        fit <- fit.hybridPLS(
          w_tr, y_tr, n_iter = n_comp, lambda = lam,
          validation_data = list(W_test = w_val, y_test = y_val)
        )
        mses[k] <- min(fit$validation_rmse^2)
      }, error = function(e) mses[k] <<- Inf)
    }
    grid_results[g] <- mean(mses)
  }
  lambda_grid[[which.min(grid_results)]]
}

one_rep <- function(rep_id) {
  set.seed(rep_id * 555)
  sim_data <- generate_strict_orthogonal_data(n_train + n_test)
  processed <- split_and_normalize.all(
    sim_data$W, sim_data$y, train_ratio = n_train / (n_train + n_test)
  )
  W_train <- processed$predictor_train
  W_test <- processed$predictor_test
  y_train <- processed$response_train
  y_test <- processed$response_test
  scale_y <- sd(y_test)
  sst <- sum((y_test - mean(y_test))^2)

  best_lam <- tune_hybrid_pls(W_train, y_train, n_folds = 5, n_comp = n_comp_max)
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
    W_train, y_train, W_test, y_test, n_comp_max
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
write.csv(ours, file.path(out_dir, "replicate_figure4_scenario1_results.csv"), row.names = FALSE)
saveRDS(ours, file.path(out_dir, "replicate_figure4_scenario1_results.rds"))

summarize_rmse <- function(df) {
  do.call(rbind, lapply(split(df, list(df$Method, df$Component), drop = TRUE), function(s) {
    data.frame(
      Method = s$Method[1],
      Component = s$Component[1],
      n = nrow(s),
      mean = mean(s$Scaled_RMSE),
      median = median(s$Scaled_RMSE),
      se = sd(s$Scaled_RMSE) / sqrt(nrow(s)),
      stringsAsFactors = FALSE
    )
  }))
}

summarize_corr <- function(df) {
  pcr <- df[df$Method == "HybridPCR", ]
  do.call(rbind, lapply(split(pcr, pcr$Component), function(s) {
    data.frame(
      Component = s$Component[1],
      mean_abs_F1F2 = mean(abs(s$Corr_F1_F2), na.rm = TRUE),
      mean_abs_F1Z = mean(abs(s$Corr_F1_Z), na.rm = TRUE),
      mean_abs_F2Z = mean(abs(s$Corr_F2_Z), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

old <- read.csv(old_csv)
cat("\n=== Scaled RMSE means (old 100 MC vs ours) ===\n")
sum_old <- summarize_rmse(old)
sum_new <- summarize_rmse(ours)
cmp <- merge(
  sum_old, sum_new,
  by = c("Method", "Component"),
  suffixes = c("_old", "_new")
)
cmp$d_mean <- cmp$mean_new - cmp$mean_old
cmp$rel <- abs(cmp$d_mean) / pmax(cmp$mean_old, 1e-12)
print(cmp[order(cmp$Component, cmp$Method),
          c("Method", "Component", "mean_old", "mean_new", "d_mean", "rel",
            "median_old", "median_new")],
      digits = 4, row.names = FALSE)

cat("\n=== PCR |cor| means (panel b) ===\n")
c_old <- summarize_corr(old)
c_new <- summarize_corr(ours)
cc <- merge(c_old, c_new, by = "Component", suffixes = c("_old", "_new"))
print(cc[order(cc$Component), ], digits = 4, row.names = FALSE)

# Key paper claims for Scenario 1 / Fig 4
pls1_old <- sum_old$mean[sum_old$Method == "HybridPLS" & sum_old$Component == 1]
pcr1_old <- sum_old$mean[sum_old$Method == "HybridPCR" & sum_old$Component == 1]
pls1_new <- sum_new$mean[sum_new$Method == "HybridPLS" & sum_new$Component == 1]
pcr1_new <- sum_new$mean[sum_new$Method == "HybridPCR" & sum_new$Component == 1]

cat("\n=== Checks ===\n")
cat(sprintf("Paper text ~ PLS 0.25 / PCR 0.66 at L=1\n"))
cat(sprintf("Old CSV:  PLS %.3f / PCR %.3f\n", pls1_old, pcr1_old))
cat(sprintf("Ours:     PLS %.3f / PCR %.3f\n", pls1_new, pcr1_new))
cat(sprintf("PLS << PCR at L=1: %s\n", pls1_new < 0.5 * pcr1_new))
# after L=2, methods should be close
pls2 <- sum_new$mean[sum_new$Method == "HybridPLS" & sum_new$Component == 2]
pcr2 <- sum_new$mean[sum_new$Method == "HybridPCR" & sum_new$Component == 2]
cat(sprintf("Converge by L=2 (both <0.15): PLS=%.3f PCR=%.3f -> %s\n",
            pls2, pcr2, pls2 < 0.15 && pcr2 < 0.15))
# first two PCR components nearly perfect cross-modality cor
c1 <- c_new[c_new$Component == 1, ]
c2 <- c_new[c_new$Component == 2, ]
high_cor <- all(c(c1$mean_abs_F1F2, c1$mean_abs_F1Z, c1$mean_abs_F2Z,
                  c2$mean_abs_F1F2, c2$mean_abs_F1Z, c2$mean_abs_F2Z) > 0.8)
cat(sprintf("PCR |cor| high for L=1,2: %s\n", high_cor))

# match old within ~15% relative on L=1 means (MC noise)
match_l1 <- abs(pls1_new - pls1_old) / pls1_old < 0.20 &&
  abs(pcr1_new - pcr1_old) / pcr1_old < 0.20
overall <- match_l1 && (pls1_new < 0.5 * pcr1_new) && high_cor
cat(sprintf("\nFIGURE 4 OVERALL vs old result.csv: %s\n",
            if (overall) "PASS" else "REVIEW"))
cat("Saved replicate_figure4_scenario1_results.{rds,csv}\n")
