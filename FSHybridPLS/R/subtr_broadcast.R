# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Broadcast Subtraction for Matrices
#'
#' Subtracts a single observation (vector) from every row of a matrix.
#'
#' @param input A numeric matrix.
#' @param other A numeric vector or a 1-row matrix.
#' @param alpha A scalar multiplier (default 1).
#' @return A matrix with the broadcasted subtraction applied.
#' @export
subtr_broadcast <- function(input, other, alpha = 1) {
  # Reuse add_broadcast with negated alpha
  add_broadcast(input, other, (-1 * alpha))
}
