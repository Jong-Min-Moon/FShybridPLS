# Figure 3 (fig:beta_sensitivity) via 6_2_beta_estimation.R logic
# Paper: n in {200,1000,3000}, lambda grid {0,0.001,0.1}^2, L=10, mean over 100 MC
# Here: full grid; MC counts chosen for runtime (paper used 100)
suppressPackageStartupMessages({
  library(FSHybridPLS)
  library(fda)
})

set.seed(20260126)
n_comp <- 10L
lam_grid <- c(0, 0.001, 0.1)
# Full paper n's; MC reduced for large n to keep runtime practical
mc_by_n <- c(`200` = 50L, `1000` = 20L, `3000` = 10L)

cat("=== Figure 3 / Section 6.2 beta sensitivity ===\n")
cat("Package:", as.character(packageVersion("FSHybridPLS")), "\n")
cat("L =", n_comp, "\n")
cat("lambda grid:", paste(lam_grid, collapse = ", "), "\n")
cat("MC by n:", paste(sprintf("n=%s:%s", names(mc_by_n), mc_by_n), collapse = "; "), "\n\n")

generate_complex_data <- function(n_sample = 200) {
  eval_points <- seq(0, 1, length.out = 100)
  basis <- create.bspline.basis(c(0, 1), 20)
  nbasis <- basis$nbasis
  beta1_fd <- Data2fd(eval_points, 2 * eval_points * sin(3 * pi * eval_points), basis)
  beta2_fd <- Data2fd(eval_points, 2 * exp(-10 * (eval_points - 0.5)^2), basis)
  beta_z_true <- c(1.5, -1.0)
  coefs1 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- matrix(rnorm(n_sample * nbasis), nbasis, n_sample)
  coefs2 <- 0.6 * coefs2 + 0.4 * coefs1
  fd1 <- fd(coefs1, basis)
  fd2 <- fd(coefs2, basis)
  Z <- t(coefs1[1:2, ]) + matrix(rnorm(n_sample * 2, sd = 0.5), n_sample, 2)
  y_clean <- as.vector(fda::inprod(fd1, beta1_fd) + fda::inprod(fd2, beta2_fd) + Z %*% beta_z_true)
  y <- y_clean + rnorm(n_sample, sd = sd(y_clean) * 0.05)
  list(
    W = predictor_hybrid(Z, list(fd1, fd2), eval_points),
    y = y,
    truth = list(beta_fd = list(beta1_fd, beta2_fd), beta_z = beta_z_true)
  )
}

one_cell <- function(n_sample, l1, l2, rep_id) {
  set.seed(rep_id + 1000L * n_sample + as.integer(l1 * 1e6) + as.integer(l2 * 1e3))
  sim <- generate_complex_data(n_sample)
  norm_t1 <- sqrt(fda::inprod(sim$truth$beta_fd[[1]], sim$truth$beta_fd[[1]]))
  norm_t2 <- sqrt(fda::inprod(sim$truth$beta_fd[[2]], sim$truth$beta_fd[[2]]))
  norm_ts <- sqrt(sum(sim$truth$beta_z^2))
  fit <- fit.hybridPLS(sim$W, sim$y, n_iter = n_comp, lambda = c(l1, l2))
  beta_est <- fit$beta[[n_comp]]
  d1 <- minus.fd(sim$truth$beta_fd[[1]], beta_est$functional_list[[1]])
  d2 <- minus.fd(sim$truth$beta_fd[[2]], beta_est$functional_list[[2]])
  ds <- sim$truth$beta_z - as.vector(beta_est$Z)
  c(
    Err_Beta1 = sqrt(fda::inprod(d1, d1)) / norm_t1,
    Err_Beta2 = sqrt(fda::inprod(d2, d2)) / norm_t2,
    Err_Scalar = sqrt(sum(ds^2)) / norm_ts
  )
}

rows <- list()
for (n_chr in names(mc_by_n)) {
  n_sample <- as.integer(n_chr)
  n_mc <- mc_by_n[[n_chr]]
  cat(sprintf("\n--- n=%d, MC=%d ---\n", n_sample, n_mc))
  for (l1 in lam_grid) {
    for (l2 in lam_grid) {
      mat <- sapply(seq_len(n_mc), function(r) one_cell(n_sample, l1, l2, r))
      mu <- rowMeans(mat)
      rows[[length(rows) + 1]] <- data.frame(
        n = n_sample, lam1 = l1, lam2 = l2, n_mc = n_mc,
        mean_Err_Beta1 = mu["Err_Beta1"],
        mean_Err_Beta2 = mu["Err_Beta2"],
        mean_Err_Scalar = mu["Err_Scalar"]
      )
      cat(sprintf("  (l1=%g,l2=%g) B1=%.3f B2=%.3f Sc=%.3f\n",
                  l1, l2, mu["Err_Beta1"], mu["Err_Beta2"], mu["Err_Scalar"]))
    }
  }
}

tab <- do.call(rbind, rows)
rownames(tab) <- NULL

# Qualitative checks from paper caption/text:
# - regularization matters (errors vary with lambda)
# - interaction between lambda1 and lambda2
# - larger n tends to reduce errors (for similar lambda)
check_varies <- function(n0, metric) {
  sub <- tab[tab$n == n0, ]
  diff(range(sub[[metric]])) > 0.02
}
check_n_helps <- function(metric) {
  # compare mean over grid at n=200 vs n=1000
  m200 <- mean(tab[tab$n == 200, metric])
  m1000 <- mean(tab[tab$n == 1000, metric])
  m1000 < m200
}

cat("\n=== Aggregated means (Figure 3 cells) ===\n")
print(tab, digits = 3, row.names = FALSE)

cat("\n=== Checks vs paper Figure 3 description ===\n")
v_b1 <- check_varies(200, "mean_Err_Beta1")
v_b2 <- check_varies(200, "mean_Err_Beta2")
v_sc <- check_varies(200, "mean_Err_Scalar")
# interaction: for fixed lam1, changing lam2 changes error
sub200 <- tab[tab$n == 200, ]
interact <- FALSE
for (l1 in lam_grid) {
  s <- sub200[sub200$lam1 == l1, ]
  if (diff(range(s$mean_Err_Beta1)) > 0.01 || diff(range(s$mean_Err_Beta2)) > 0.01) interact <- TRUE
}
n_helps <- check_n_helps("mean_Err_Beta1") && check_n_helps("mean_Err_Beta2")

cat(sprintf("Errors vary with lambda (n=200): B1=%s B2=%s Sc=%s\n", v_b1, v_b2, v_sc))
cat(sprintf("Lambda1/Lambda2 interaction visible: %s\n", interact))
cat(sprintf("Larger n reduces avg error (200->1000): %s\n", n_helps))
cat(sprintf("\nFIGURE 3 OVERALL: %s\n",
            if (v_b1 && v_b2 && v_sc && interact && n_helps) "PASS (pattern match)" else "REVIEW"))

out_dir <- "C:/Users/jongmin/Downloads/FShybridPLS-main/FShybridPLS-main/numerical_studies"
saveRDS(list(table = tab, n_comp = n_comp, mc_by_n = mc_by_n),
        file.path(out_dir, "replicate_figure3_beta_results.rds"))
# also write CSV for easy viewing
write.csv(tab, file.path(out_dir, "replicate_figure3_beta_results.csv"), row.names = FALSE)
cat("Saved results to numerical_studies/replicate_figure3_beta_results.{rds,csv}\n")
