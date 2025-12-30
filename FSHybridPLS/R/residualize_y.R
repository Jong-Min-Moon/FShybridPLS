# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Residualize the response vector
#'
#' Subtracts the projection of y onto rho from y.
#'
#' @param y The current response vector.
#' @param rho The current PLS score vector.
#' @param nu The regression coefficient.
#'
#' @return The residualized response vector.
#' @export
residualize_y <- function(y, rho, nu) {
  y_next <- y - nu * rho
  return(y_next)
}
