# Replicate Figure 2 (Y vs PLS score correlations) under refactored FSHybridPLS
# Paper: Mixed (lambda=0.1,10), mean |cor(Y, rho^[l])| for l=1..5 over 100 MC, error bars = 1 SD
suppressPackageStartupMessages({
  library(FSHybridPLS)
  library(fda)
})

n_mc <- 100L
n_comp <- 5L
lam <- c(0.1, 10.0)  # Mixed
set.seed(20260126)

cat("=== Figure 2 replication ===\n")
cat("Package:", as.character(packageVersion("FSHybridPLS")), "\n")
cat("n_mc =", n_mc, ", L =", n_comp, ", lambda =", paste(lam, collapse = ","), "\n\n")

generate_torture_data <- function(n_sample = 100) {
  t1 <- rnorm(n_sample, sd = 10)
  t2 <- rnorm(n_sample, sd = 0.1)
  eval_points <- seq(0, 1, length.out = 100)
  basis <- create.bspline.basis(c(0, 1), 15)
  fd1 <- Data2fd(eval_points, t(outer(t1, sin(2 * pi * eval_points))), basis)
  noise_func <- matrix(rnorm(n_sample * 100, sd = 0.01), 100, n_sample)
  fd2 <- Data2fd(eval_points, t(outer(t2, sin(10 * pi * eval_points)) + noise_func), basis)
  P <- 50
  A <- matrix(runif(P * 2, min = -1, max = 1), nrow = P, ncol = 2)
  Z <- cbind(t1, t2) %*% t(A) + matrix(rnorm(n_sample * P, sd = 1), n_sample, P)
  y <- 0.5 * t1 + 10 * t2 + rnorm(n_sample, sd = 1)
  list(W = predictor_hybrid(Z, list(fd1, fd2), eval_points), y = y)
}

cors_mat <- matrix(NA_real_, n_mc, n_comp)
for (r in seq_len(n_mc)) {
  data <- generate_torture_data(100)
  processed <- split_and_normalize.all(data$W, data$y, train_ratio = 0.99)
  W <- processed$predictor_train
  y <- processed$response_train
  fit <- fit.hybridPLS(W, y, n_iter = n_comp, lambda = lam)
  Rho <- do.call(cbind, fit$rho)
  cors_mat[r, ] <- abs(as.numeric(cor(y, Rho)))
  if (r %% 20 == 0) cat("rep", r, "\n")
}

mu <- colMeans(cors_mat)
sdev <- apply(cors_mat, 2, sd)

res <- data.frame(
  component = seq_len(n_comp),
  mean_abs_cor = mu,
  sd = sdev
)
print(res, digits = 3, row.names = FALSE)

# Paper qualitative claims for Figure 2
decaying <- all(diff(mu) < 0)
first_largest <- mu[1] == max(mu)
second_positive <- mu[2] > 0.05

cat("\n=== Checks vs paper Figure 2 description ===\n")
cat(sprintf("Component 1 has highest |cor|: %s (%.3f)\n", first_largest, mu[1]))
cat(sprintf("Monotone decay l=1..5: %s\n", decaying))
cat(sprintf("Component 2 still non-negligible: %s (%.3f)\n", second_positive, mu[2]))
cat(sprintf("\nFIGURE 2 OVERALL: %s\n",
            if (first_largest && decaying && second_positive) "PASS (pattern match)" else "REVIEW"))

# Save for optional plotting
out <- list(summary = res, all = cors_mat, lambda = lam, n_mc = n_mc)
saveRDS(out, file.path(tempdir(), "fig2_y_score_cors.rds"))
cat("Saved:", file.path(tempdir(), "fig2_y_score_cors.rds"), "\n")
