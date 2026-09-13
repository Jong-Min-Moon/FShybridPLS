#' Compute scalar regression coefficient (nu)
#'
#' @param y Response vector.
#' @param rho PLS score vector.
#' @return Scalar regression coefficient.
#' @keywords internal
get_nu <- function(y, rho) {
  denom <- sum(rho * rho)
  if (!is.finite(denom) || denom <= .Machine$double.eps) {
    stop("Score vector 'rho' has zero norm; cannot compute nu.", call. = FALSE)
  }
  sum(y * rho) / denom
}
