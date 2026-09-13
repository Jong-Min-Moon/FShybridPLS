# Backward-compatible aliases for paper replication scripts.
# These preserve the original published function names without changing
# the numerical algorithm (wrappers around the refactored implementations).
# Note: do not export split.all — it clashes with base::split S3 generic.

#' @rdname fit_hybridPLS
#' @export
fit.hybridPLS <- function(W, y, n_iter, lambda, validation_data = NULL) {
  fit_hybridPLS(W, y, n_iter, lambda, validation_data)
}

#' @rdname inprod_predictor_hybrid
#' @export
inprod.predictor_hybrid <- function(xi_1, xi_2 = NULL) {
  inprod_predictor_hybrid(xi_1, xi_2)
}

#' @rdname inprod_pen_predictor_hybrid
#' @export
inprod_pen.predictor_hybrid <- function(xi_1, xi_2 = NULL, lambda) {
  inprod_pen_predictor_hybrid(xi_1, xi_2, lambda)
}

#' @rdname split_and_normalize_all
#' @export
split_and_normalize.all <- function(W_hybrid, response, train_ratio) {
  split_and_normalize_all(W_hybrid, response, train_ratio)
}
