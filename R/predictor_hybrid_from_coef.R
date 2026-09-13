#' Construct a Single-Sample Predictor Hybrid Object from Coefficients
#'
#' @param format A \code{predictor_hybrid} template providing structure/basis info.
#' @param coef A numeric vector of functional then scalar coefficients.
#'
#' @return A single-sample \code{predictor_hybrid} object.
#' @keywords internal
#' @importFrom fda fd
predictor_hybrid_from_coef <- function(format, coef) {
  assert_predictor_hybrid(format, "format")
  coef <- as.numeric(coef)

  basis_counts <- format$n_basis_list
  if (length(unique(basis_counts)) != 1) {
    stop("All functional predictors must have the same number of basis functions (M).", call. = FALSE)
  }

  M <- basis_counts[1]
  K <- format$n_functional
  expected_len <- K * M + format$n_scalar
  if (length(coef) != expected_len) {
    stop(sprintf("'coef' must have length %d.", expected_len), call. = FALSE)
  }

  new_fun <- vector("list", K)
  for (ii in seq_len(K)) {
    start_idx <- (ii - 1) * M + 1
    end_idx <- ii * M
    new_fun[[ii]] <- fd(
      coef = as.matrix(coef[start_idx:end_idx]),
      basisobj = format$functional_list[[ii]]$basis
    )
  }

  new_Z <- matrix(coef[(K * M + 1):length(coef)], nrow = 1)
  if (!is.null(colnames(format$Z))) {
    colnames(new_Z) <- colnames(format$Z)
  }

  predictor_hybrid(
    Z = new_Z,
    functional_list = new_fun,
    eval_point = format$eval_point
  )
}
