# Generated from create-rhello.Rmd: do not edit by hand

#' Multiply a predictor_hybrid object by a scalar
#'
#' Performs scalar multiplication of both scalar and functional components 
#' of a `predictor_hybrid` object. Functional predictors are scaled using 
#' `times.fd()` from the `fda` package.
#'
#' @param input A `predictor_hybrid` object.
#' @param scalar A numeric value to multiply all components by.
#'
#' @return A new `predictor_hybrid` object scaled by `scalar`.
#' @export
scalar_mul.predictor_hybrid <- function(input, scalar) {
  new_functional_list <- lapply(
    input$functional_list,
    function(fd_obj) times.fd(scalar, fd_obj)
  )
  
  new_Z <- scalar * input$Z
  
  predictor_hybrid(
    Z = new_Z,
    functional_list = new_functional_list,
    jacobian_list = input$jacobian_list,
    n_basis_list = input$n_basis_list,
    n_sample = input$n_sample,
    n_functional = input$n_functional,
    n_scalar = input$n_scalar
  )
}
