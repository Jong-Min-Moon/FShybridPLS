# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Multiply a predictor_hybrid object by a scalar
#'
#' Performs scalar multiplication on both the scalar and functional components
#' of a `predictor_hybrid` object. Functional predictors are scaled using
#' `times.fd()` from the `fda` package.
#'
#' @param input A `predictor_hybrid` object.
#' @param scalar A numeric value to multiply all components by.
#'
#' @return A new `predictor_hybrid` object scaled by `scalar`.
#' @export
scalar_mul.predictor_hybrid <- function(input, scalar) {
  if (!inherits(input, "predictor_hybrid")) {
    stop("Input must be of class 'predictor_hybrid'.")
  }
  if (!is.numeric(scalar) || length(scalar) != 1) {
    stop("Scalar must be a single numeric value.")
  }

  # Scale functional components
  new_functional_list <- lapply(input$functional_list, function(fd_obj) {
    times.fd(scalar, fd_obj)
  })

  # Scale scalar predictors
  new_Z <- scalar * input$Z

  # Reconstruct hybrid object
  # Passing eval_point from the original input to maintain metadata
  predictor_hybrid(
    Z = new_Z,
    functional_list = new_functional_list,
    eval_point = input$eval_point
  )
}
