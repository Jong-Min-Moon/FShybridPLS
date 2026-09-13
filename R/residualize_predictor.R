#' Residualize the hybrid predictor
#'
#' Updates each observation in W by subtracting the projection onto rho
#' (`rho_i * delta`) using vectorized linear algebra.
#'
#' @param W The current `predictor_hybrid` object.
#' @param rho The current PLS score vector.
#' @param delta The hybrid regression coefficient (single-sample object).
#'
#' @return The residualized `predictor_hybrid` object.
#' @keywords internal
#' @importFrom fda fd
residualize_predictor <- function(W, rho, delta) {
  rho <- as.vector(rho)
  if (length(rho) != W$n_sample) {
    stop("'rho' length must equal W$n_sample.", call. = FALSE)
  }
  if (!inherits(delta, "predictor_hybrid") || delta$n_sample != 1) {
    stop("'delta' must be a single-sample predictor_hybrid.", call. = FALSE)
  }

  new_Z <- W$Z - rho %*% delta$Z
  new_fun <- lapply(seq_len(W$n_functional), function(k) {
    old <- W$functional_list[[k]]
    new_coefs <- old$coefs - (delta$functional_list[[k]]$coefs %*% t(rho))
    fd(coef = new_coefs, basisobj = old$basis, fdnames = old$fdnames)
  })

  out <- W
  out$Z <- new_Z
  out$functional_list <- new_fun
  out
}
