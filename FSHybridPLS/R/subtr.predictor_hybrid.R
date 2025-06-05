# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Subtract two predictor_hybrid objects
#'
#' Performs element-wise subtraction of two `predictor_hybrid` objects.
#' Internally uses `add.predictor_hybrid(input, other, alpha = -1)`.
#'
#' @param input A `predictor_hybrid` object.
#' @param other Another `predictor_hybrid` object to subtract.
#' @param alpha A scalar multiplier applied to `other` before subtraction (default is 1).
#'
#' @return A new `predictor_hybrid` object representing the result of subtraction.
#' @export
subtr.predictor_hybrid <- function(input, other, alpha = 1) {
  add.predictor_hybrid(input, other, alpha = -alpha)
}

