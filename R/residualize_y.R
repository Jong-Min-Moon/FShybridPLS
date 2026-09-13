#' Residualize the response vector
#'
#' @param y Response vector.
#' @param rho Score vector.
#' @param nu Regression coefficient.
#' @return Residualized response.
#' @keywords internal
residualize_y <- function(y, rho, nu) {
  as.numeric(y - nu * rho)
}
