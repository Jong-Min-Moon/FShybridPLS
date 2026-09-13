# Build frame data for the Hybrid PLS component-extraction animation.
# Richer predictive functional shapes; original predictor–Y correlations included.

suppressPackageStartupMessages({
  library(FSHybridPLS)
  library(fda)
})

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org")
}
library(jsonlite)

set.seed(7)
out_dir <- "D:/GitHub/jongmin/FShybridPLS/animation"
n <- 150L
n_show <- 8L
n_comp <- 2L
eval_points <- seq(0, 1, length.out = 90)
basis <- create.bspline.basis(c(0, 1), 18)
t <- eval_points

# Predictive coefficient shapes (intentionally different / non-flat)
beta1 <- exp(-((t - 0.28) / 0.12)^2) - 0.45 * exp(-((t - 0.72) / 0.14)^2)
# Delayed multiphasic predictive shape for X2 (visually distinct from X1)
beta2 <- exp(-((t - 0.55) / 0.10)^2) - 0.7 * exp(-((t - 0.22) / 0.11)^2) +
  0.55 * sin(4 * pi * t) * exp(-((t - 0.65) / 0.35)^2)
beta1 <- as.numeric(scale(beta1))
beta2 <- as.numeric(scale(beta2))

# High-variance nuisance shapes (different geometry)
noise1 <- sin(6 * pi * t) + 0.4 * cos(2 * pi * t)
noise2 <- cos(5 * pi * t) - 0.5 * sin(3 * pi * t)
noise1 <- as.numeric(scale(noise1))
noise2 <- as.numeric(scale(noise2))

u_signal <- rnorm(n, sd = 1)
y_clean <- 2.2 * u_signal
y_raw <- y_clean + rnorm(n, sd = sd(y_clean) * 0.1)

u_noise_raw <- rnorm(n, sd = 1)
u_noise <- as.numeric(scale(residuals(lm(u_noise_raw ~ y_raw)))) * 2.6

X1 <- outer(u_signal, beta1) + outer(u_noise, noise1) +
  matrix(rnorm(n * length(t), sd = 0.15), n)
# Stronger predictive loading on X2 so ξ₂ is not tiny relative to nuisance
X2 <- outer(1.35 * u_signal, beta2) + outer(0.85 * u_noise, noise2) +
  matrix(rnorm(n * length(t), sd = 0.15), n)

fd1 <- Data2fd(t, t(X1), basis)
fd2 <- Data2fd(t, t(X2), basis)

Z <- cbind(
  u_noise + rnorm(n, 0.12),
  u_noise + rnorm(n, 0.12),
  0.55 * u_noise + 0.45 * u_signal + rnorm(n, 0.12),
  u_signal + rnorm(n, 0.12)
)
colnames(Z) <- c("Z1", "Z2", "Z3", "Z4")

W <- predictor_hybrid(Z, list(fd1, fd2), t)
y <- as.numeric(scale(y_raw))

# Mild penalty so functional weights keep shape detail
fit <- fit_hybridPLS(W, y, n_iter = n_comp, lambda = c(0.02, 0.02))

# ---- Original predictor–response associations ----
cor_fun <- function(Xmat, yvec) {
  apply(Xmat, 2, function(col) {
    if (sd(col) < 1e-10) return(0)
    as.numeric(cor(col, yvec))
  })
}
orig_cor_x1 <- cor_fun(X1, y)
orig_cor_x2 <- cor_fun(X2, y)
orig_cor_z <- as.numeric(cor(Z, y))
names(orig_cor_z) <- colnames(Z)

# Summary: how much association lives in original predictors vs one score
# Use mean absolute correlation for functionals + abs scalar cors
orig_assoc_summary <- c(
  mean_abs_cor_X1 = mean(abs(orig_cor_x1)),
  mean_abs_cor_X2 = mean(abs(orig_cor_x2)),
  mean_abs_cor_Z = mean(abs(orig_cor_z)),
  max_abs_cor_Z = max(abs(orig_cor_z))
)

eval_curves <- function(W_obj, ids = NULL) {
  if (is.null(ids)) ids <- seq_len(W_obj$n_sample)
  out <- list()
  for (k in seq_along(W_obj$functional_list)) {
    mat <- eval.fd(W_obj$eval_point, W_obj$functional_list[[k]][ids])
    out[[k]] <- lapply(seq_along(ids), function(j) as.numeric(mat[, j]))
  }
  out
}

eval_xi <- function(xi_obj) {
  list(
    functional = lapply(xi_obj$functional_list, function(fd_obj) {
      as.numeric(eval.fd(xi_obj$eval_point, fd_obj))
    }),
    scalar = as.numeric(xi_obj$Z[1, ])
  )
}

# Subjects ordered by Y for an explicit response legend
ord <- order(y)
show_ids <- ord[as.integer(round(seq(1, n, length.out = n_show)))]

y_show <- y[show_ids]
y01 <- (y_show - min(y)) / (max(y) - min(y) + 1e-8)

curves0 <- eval_curves(fit$W[[1]], show_ids)
curves1 <- eval_curves(fit$W[[2]], show_ids)
curves2 <- eval_curves(fit$W[[3]], show_ids)
Z0 <- fit$W[[1]]$Z[show_ids, , drop = FALSE]
Z1 <- fit$W[[2]]$Z[show_ids, , drop = FALSE]
Z2 <- fit$W[[3]]$Z[show_ids, , drop = FALSE]

y_res <- vector("list", n_comp + 1)
y_res[[1]] <- y
y_now <- y
for (l in seq_len(n_comp)) {
  y_now <- y_now - fit$nu[[l]] * fit$rho[[l]]
  y_res[[l + 1]] <- y_now
}

cors <- vapply(seq_len(n_comp), function(l) {
  abs(cor(fit$rho[[l]], y_res[[l]]))
}, numeric(1))
# Also absolute correlation of first score with original Y
cors_vs_orig_y <- vapply(seq_len(n_comp), function(l) {
  abs(cor(fit$rho[[l]], y))
}, numeric(1))

# Sign-align xi so it matches the sign of original correlations (visual overlay)
xi1 <- eval_xi(fit$xi[[1]])
# Flip if needed so xi1 resembles orig_cor_x1
if (cor(xi1$functional[[1]], orig_cor_x1) < 0) {
  xi1$functional <- lapply(xi1$functional, `-`)
  xi1$scalar <- -xi1$scalar
  # rho sign flips with xi; keep payload consistent by flipping stored rho for display
  flip1 <- -1
} else {
  flip1 <- 1
}

payload <- list(
  meta = list(
    title = "Hybrid PLS: supervised basis extraction",
    tagline = "One joint direction (X1 + X2 + scalars) collapses into one score that keeps the predictor–Y association",
    n = n,
    n_show = length(show_ids),
    n_comp = n_comp,
    t = as.numeric(t),
    scalar_names = colnames(Z),
    cors = as.numeric(cors),
    cors_vs_orig_y = as.numeric(cors_vs_orig_y),
    orig_assoc_summary = as.list(orig_assoc_summary),
    note = "Synthetic data with distinct predictive shapes in X1/X2 plus high-variance nuisance orthogonal to Y"
  ),
  original_association = list(
    cor_X1 = as.numeric(orig_cor_x1),
    cor_X2 = as.numeric(orig_cor_x2),
    cor_Z = as.numeric(orig_cor_z),
    true_beta1 = as.numeric(beta1),
    true_beta2 = as.numeric(beta2)
  ),
  subjects = lapply(seq_along(show_ids), function(j) {
    list(
      id = as.integer(show_ids[j]),
      y = y_show[j],
      y01 = y01[j],
      stages = list(
        list(curves = list(curves0[[1]][[j]], curves0[[2]][[j]]), scalars = as.numeric(Z0[j, ])),
        list(curves = list(curves1[[1]][[j]], curves1[[2]][[j]]), scalars = as.numeric(Z1[j, ])),
        list(curves = list(curves2[[1]][[j]], curves2[[2]][[j]]), scalars = as.numeric(Z2[j, ]))
      )
    )
  }),
  components = lapply(seq_len(n_comp), function(l) {
    xi <- eval_xi(fit$xi[[l]])
    rho <- as.numeric(fit$rho[[l]])
    if (l == 1 && flip1 < 0) {
      xi$functional <- lapply(xi$functional, `-`)
      xi$scalar <- -xi$scalar
      rho <- -rho
    }
    list(
      l = l,
      xi_functional = xi$functional,
      xi_scalar = xi$scalar,
      rho = rho,
      y_before = as.numeric(y_res[[l]]),
      y_after = as.numeric(y_res[[l + 1]]),
      cor_abs = cors[l],
      cor_vs_orig_y = cors_vs_orig_y[l],
      nu = fit$nu[[l]]
    )
  }),
  scatter_idx = as.integer(seq(1, n, length.out = 70))
)

json_path <- file.path(out_dir, "animation_data.json")
js_path <- file.path(out_dir, "animation_data.js")
write_json(payload, json_path, auto_unbox = TRUE, pretty = TRUE, digits = 6)
json_txt <- as.character(toJSON(payload, auto_unbox = TRUE, digits = 6))
writeLines(paste0("window.ANIM_DATA = ", json_txt, ";"), js_path)

cat("Wrote", json_path, "\n")
cat("Wrote", js_path, "\n")
cat("Orig mean |cor| X1/X2/Z:",
    round(orig_assoc_summary, 3), "\n")
cat("|cor(rho,y)| by component:", paste(round(cors_vs_orig_y, 3), collapse = ", "), "\n")
cat("Open:", file.path(out_dir, "index.html"), "\n")
