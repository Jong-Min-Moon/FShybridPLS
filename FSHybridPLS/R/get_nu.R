# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute scalar regression coefficient (nu)
#'
#' Calculates the projection coefficient of the response vector y onto the score vector rho.
#'
#' @param y A numeric vector of response values.
#' @param rho A numeric vector of PLS scores.
#'
#' @return A scalar regression coefficient.
#' @export
get_nu <- function(y, rho) {
  # Compute numerator: inner product of y and rho
  # Compute denominator: squared norm of rho
  nu <- sum(y * rho) / sum(rho * rho)
  return(nu)
}
