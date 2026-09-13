#' Compute hybrid regression coefficient (delta)
#'
#' Weighted average of predictor observations with weights `rho`.
#'
#' @param W A `predictor_hybrid` object.
#' @param rho A numeric vector of PLS scores.
#' @return A single-sample `predictor_hybrid` object.
#' @keywords internal
#' @importFrom fda fd
get_delta <- function(W, rho) {
  rho <- as.vector(rho)
  denom <- sum(rho * rho)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    stop("Score vector 'rho' has zero norm; cannot compute delta.", call. = FALSE)
  }

  new_Z <- matrix(as.numeric(t(rho) %*% W$Z) / denom, nrow = 1)
  if (!is.null(colnames(W$Z))) colnames(new_Z) <- colnames(W$Z)

  new_fun <- lapply(seq_len(W$n_functional), function(k) {
    old <- W$functional_list[[k]]
    # coefs: M x n; weighted mean column = (coefs %*% rho) / denom
    new_coefs <- as.matrix((old$coefs %*% rho) / denom)
    fd(coef = new_coefs, basisobj = old$basis, fdnames = old$fdnames)
  })

  predictor_hybrid(
    Z = new_Z,
    functional_list = new_fun,
    eval_point = W$eval_point
  )
}
