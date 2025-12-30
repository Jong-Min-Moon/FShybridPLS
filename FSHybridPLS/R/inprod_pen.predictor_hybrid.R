# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Penalized Inner product between two predictor_hybrid objects
#'
#' Computes the inner product between two `predictor_hybrid` objects,
#' adding a roughness penalty term based on the second derivative.
#'
#' @param xi_1 A `predictor_hybrid` object.
#' @param xi_2 Another `predictor_hybrid` object. If missing, defaults to `xi_1`.
#' @param lambda A numeric vector of smoothing parameters (one for each functional predictor).
#'
#' @return A numeric vector of inner products.
#' @export
inprod_pen.predictor_hybrid <- function(xi_1, xi_2 = NULL, lambda) {
  
  # Compute the standard inner product first
  base_inprod <- inprod.predictor_hybrid(xi_1, xi_2)
  
  # Handle self-reference for xi_2
  if (is.null(xi_2)) xi_2 <- xi_1
  
  f1 <- xi_1$functional_list
  f2 <- xi_2$functional_list
  
  # Add penalty terms for each functional predictor
  # The penalty is lambda * <D^2 f1, D^2 f2>
  for (j in seq_len(xi_1$n_functional)) {
    base_inprod <- base_inprod + lambda[j] * fda::inprod(
      fdobj1 = f1[[j]],
      fdobj2 = f2[[j]],
      Lfdobj1 = 2, # Second derivative
      Lfdobj2 = 2  # Second derivative
    )
  }
  
  return(as.numeric(base_inprod))
}
