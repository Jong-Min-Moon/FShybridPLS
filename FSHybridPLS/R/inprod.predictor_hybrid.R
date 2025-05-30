# Generated from create-rhello.Rmd: do not edit by hand

#' Inner product between two predictor_hybrid objects (with broadcasting)
#'
#' Computes the inner product between two `predictor_hybrid` objects,
#' including both functional and scalar components. Supports broadcasting
#' when one of the inputs has a single observation. If only one input is given,
#' the inner product is computed with itself.
#'
#' @param xi_1 A `predictor_hybrid` object.
#' @param xi_2 Another `predictor_hybrid` object. If missing, defaults to `xi_1`.
#'
#' @return A numeric matrix of inner products (or a scalar if both inputs have one row).
#' @export
inprod.predictor_hybrid <- function(xi_1, xi_2 = NULL) {
  if (is.null(xi_2)) xi_2 <- xi_1
  
  # Optional: Check for basis compatibility
  # if (!is_same_basis(xi_1, xi_2)) stop("Functional predictors must share the same basis.")

  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample

  # Compute sum of pairwise functional inner products
  inprod_functional <- Reduce(`+`,
    Map(function(f1, f2) fda::inprod(f1, f2),
        xi_1$functional_list, xi_2$functional_list)
  )

  # Broadcast-aware scalar component
  Z1 <- if (n1 == 1) matrix(xi_1$Z, nrow = 1) else xi_1$Z
  Z2 <- if (n2 == 1) matrix(xi_2$Z, nrow = 1) else xi_2$Z
  inprod_scalar <- Z1 %*% t(Z2)

  result <- inprod_functional + inprod_scalar

  if (n1 == 1 && n2 == 1) as.numeric(result) else result
}

