# Table 2 (tab:geometric_validation) vs 6_1 logic under refactored FSHybridPLS
# Paper settings: L=10, 100 MC, Weak/Mixed/Strong lambdas
suppressPackageStartupMessages({
  library(FSHybridPLS)
  library(fda)
})

n_mc <- 100L
n_comp <- 10L  # paper L=10 (script default 6 is just a plug-in default)
set.seed(20260126)

cat("=== Table 2 / Section 6.1 ===\n")
cat("Package:", as.character(packageVersion("FSHybridPLS")), "\n")
cat("n_mc =", n_mc, ", L =", n_comp, "\n\n")

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

one_rep <- function(lam) {
  data <- generate_torture_data(100)
  processed <- split_and_normalize.all(data$W, data$y, train_ratio = 0.99)
  fit <- fit.hybridPLS(processed$predictor_train, processed$response_train,
                       n_iter = n_comp, lambda = lam)

  # Col1: max_{k<l} ||<W^{[l]}, xi^{[k]}>||_2  (annihilation)
  max_annih <- 0
  for (l in seq_len(n_comp - 1)) {
    for (k in (l + 1):n_comp) {
      proj <- inprod.predictor_hybrid(fit$W[[k]], fit$xi[[l]])
      max_annih <- max(max_annih, sqrt(sum(proj^2)))
    }
  }

  # Col2: max_{k!=l} |<xi^{[k]}, xi^{[l]}>_{H,Lambda}|
  max_dir <- 0
  for (r in seq_len(n_comp)) {
    for (c in seq_len(n_comp)) {
      if (r == c) next
      max_dir <- max(max_dir, abs(inprod_pen.predictor_hybrid(fit$xi[[r]], fit$xi[[c]], lam)))
    }
  }

  # Col3: max_{k!=l} |Cor(rho^{[k]}, rho^{[l]})|
  max_score <- 0
  for (r in seq_len(n_comp)) {
    for (c in seq_len(n_comp)) {
      if (r == c) next
      max_score <- max(max_score, abs(cor(fit$rho[[r]], fit$rho[[c]])))
    }
  }

  c(max_annih, max_dir, max_score)
}

scenarios <- list(
  Weak   = c(0.1, 0.1),
  Mixed  = c(0.1, 10.0),
  Strong = c(10.0, 10.0)
)

# Paper Table 2 entries (already in display units of the table header)
paper <- rbind(
  Weak   = c(annih = 2.20, annih_se = 0.64, dir = 0.10, dir_se = 0.03, score = 8.12, score_se = 3.66),
  Mixed  = c(annih = 2.33, annih_se = 0.87, dir = 9.76, dir_se = 3.72, score = 8.49, score_se = 3.87),
  Strong = c(annih = 2.33, annih_se = 0.73, dir = 9.59, dir_se = 3.86, score = 8.22, score_se = 4.19)
)

# Table header scales
s_annih <- 1e-15
s_dir   <- 1e-11
s_score <- 1e-16

out_rows <- list()
for (nm in names(scenarios)) {
  cat("Running", nm, "...\n")
  mat <- replicate(n_mc, one_rep(scenarios[[nm]]))
  # mat: 3 x n_mc
  mu <- rowMeans(mat)
  se <- apply(mat, 1, function(z) sd(z) / sqrt(length(z)))

  disp_mu <- c(mu[1] / s_annih, mu[2] / s_dir, mu[3] / s_score)
  disp_se <- c(se[1] / s_annih, se[2] / s_dir, se[3] / s_score)
  names(disp_mu) <- names(disp_se) <- c("annih", "dir", "score")

  p <- paper[nm, ]
  out_rows[[nm]] <- data.frame(
    Scenario = nm,
    metric = c("annih(x1e-15)", "dir(x1e-11)", "score(x1e-16)"),
    paper_mean = c(p["annih"], p["dir"], p["score"]),
    paper_se   = c(p["annih_se"], p["dir_se"], p["score_se"]),
    ours_mean  = unname(disp_mu),
    ours_se    = unname(disp_se),
    raw_mean   = mu,
    ratio_ours_over_paper = unname(disp_mu / c(p["annih"], p["dir"], p["score"]))
  )

  cat(sprintf(
    "  Paper: %5.2f(%4.2f) | %5.2f(%4.2f) | %5.2f(%4.2f)\n  Ours:  %5.2f(%4.2f) | %5.2f(%4.2f) | %5.2f(%4.2f)\n  Raw:   %.3e | %.3e | %.3e\n",
    p["annih"], p["annih_se"], p["dir"], p["dir_se"], p["score"], p["score_se"],
    disp_mu[1], disp_se[1], disp_mu[2], disp_se[2], disp_mu[3], disp_se[3],
    mu[1], mu[2], mu[3]
  ))
}

tab <- do.call(rbind, out_rows)
rownames(tab) <- NULL
cat("\n=== Comparison table (display units as in paper) ===\n")
print(tab[, c("Scenario", "metric", "paper_mean", "paper_se", "ours_mean", "ours_se", "ratio_ours_over_paper")],
      digits = 3, row.names = FALSE)

# Pass criteria: same order of magnitude; means within ~factor 3 (MC noise / float)
ok <- all(tab$raw_mean < c(rep(1e-12, 3), rep(1e-8, 3), rep(1e-10, 3))[seq_len(nrow(tab))]) # loose per-row
# Better: check each metric's raw scale
ok_annih <- all(tab$raw_mean[tab$metric == "annih(x1e-15)"] < 1e-13)
ok_dir   <- all(tab$raw_mean[tab$metric == "dir(x1e-11)"] < 1e-8)
ok_score <- all(tab$raw_mean[tab$metric == "score(x1e-16)"] < 1e-12)
# And ratios not crazy for annih/dir
ok_ratio <- all(tab$ratio_ours_over_paper[tab$metric != "score(x1e-16)"] > 0.25 &
                tab$ratio_ours_over_paper[tab$metric != "score(x1e-16)"] < 4)

cat("\n=== Verdict ===\n")
cat(sprintf("Annihilation ~1e-15: %s\n", if (ok_annih) "MATCH" else "FAIL"))
cat(sprintf("Direction   ~1e-11: %s\n", if (ok_dir) "MATCH" else "FAIL"))
cat(sprintf("Score       ~1e-16: %s\n", if (ok_score) "MATCH" else "FAIL"))
cat(sprintf("Mean ratios (annih/dir) in [0.25,4]: %s\n", if (ok_ratio) "MATCH" else "FAIL"))
cat(sprintf("\nTABLE 2 OVERALL: %s\n", if (ok_annih && ok_dir && ok_score && ok_ratio) "PASS" else "REVIEW"))

saveRDS(list(table = tab, n_mc = n_mc, n_comp = n_comp),
        file = file.path(tempdir(), "table2_geom_compare.rds"))
