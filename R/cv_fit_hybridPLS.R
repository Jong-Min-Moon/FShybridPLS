#' Cross-Validated Hybrid Penalized PLS
#'
#' Fits [fit_hybridPLS()] on K training folds and records mean validation RMSE
#' for each number of components. Useful for choosing `n_iter`.
#'
#' @param W A `predictor_hybrid` object.
#' @param y Numeric response vector.
#' @param n_iter Maximum number of components to consider.
#' @param lambda Numeric smoothing-parameter vector (length `W$n_functional`).
#' @param n_fold Number of CV folds (default 5).
#' @param seed Optional seed for fold assignment.
#'
#' @return A list with:
#' \describe{
#'   \item{rmse_by_component}{Mean validation RMSE for components `1:n_iter`.}
#'   \item{best_n_iter}{Component count with smallest mean RMSE.}
#'   \item{fold_rmse}{Matrix of fold-wise RMSE (`n_fold` x `n_iter`).}
#' }
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 40, n_basis = 5)
#' cv <- cv_fit_hybridPLS(sim$W, sim$y, n_iter = 3, lambda = 1e-3, n_fold = 3)
#' cv$best_n_iter
#'
#' @export
cv_fit_hybridPLS <- function(W, y, n_iter, lambda, n_fold = 5, seed = NULL) {
  assert_predictor_hybrid(W, "W")
  if (!is.numeric(y) || length(y) != W$n_sample) {
    stop("'y' must be numeric of length W$n_sample.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  folds <- create_idx_kfold(W$n_sample, n_fold)
  fold_rmse <- matrix(NA_real_, nrow = n_fold, ncol = n_iter)

  for (f in seq_along(folds)) {
    idx_train <- folds[[f]]$idx_train
    idx_valid <- folds[[f]]$idx_valid

    W_train <- subset_by_index(W, idx_train)
    W_valid <- subset_by_index(W, idx_valid)
    y_train <- y[idx_train]
    y_valid <- y[idx_valid]

    fit <- fit_hybridPLS(
      W_train,
      y_train,
      n_iter = n_iter,
      lambda = lambda,
      validation_data = list(W_test = W_valid, y_test = y_valid)
    )
    fold_rmse[f, ] <- fit$validation_rmse
  }

  rmse_by_component <- colMeans(fold_rmse)
  best_n_iter <- which.min(rmse_by_component)

  list(
    rmse_by_component = rmse_by_component,
    best_n_iter = as.integer(best_n_iter),
    fold_rmse = fold_rmse
  )
}

# Internal: subset predictor_hybrid by sample indices
subset_by_index <- function(W, idx) {
  fun_list <- lapply(W$functional_list, function(fd_obj) fd_obj[idx])
  predictor_hybrid(
    Z = W$Z[idx, , drop = FALSE],
    functional_list = fun_list,
    eval_point = W$eval_point
  )
}