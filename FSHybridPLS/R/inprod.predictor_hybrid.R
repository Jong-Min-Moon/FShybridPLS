# Generated from create-rhello.Rmd: do not edit by hand

#' Inner product between two predictor_hybrid objects (with broadcasting)
#'
#' Computes the inner product between two `predictor_hybrid` objects,
#' including both functional and scalar components. Supports broadcasting
#' when one of the inputs has a single observation.
#'
#' @param xi_1 A `predictor_hybrid` object.
#' @param xi_2 Another `predictor_hybrid` object. If missing, defaults to `xi_1`.
#'
#' @return A numeric vector of inner products, or a scalar if both inputs contain a single observation.
#' @export
inprod.predictor_hybrid <- function(xi_1, xi_2 = NULL) {
  # Handle self-inner product
  if (is.null(xi_2)) xi_2 <- xi_1 
  
  # Type checks
  if (!inherits(xi_1, "predictor_hybrid") || !inherits(xi_2, "predictor_hybrid")) {
    stop("Both inputs must be of class 'predictor_hybrid'.")
  }

  # Structural checks
  if (xi_1$n_functional != xi_2$n_functional) {
    stop("Mismatch in number of functional predictors.")
  }
  if (xi_1$n_scalar != xi_2$n_scalar) {
    stop("Mismatch in number of scalar predictors.")
  }

  # Swap so that broadcasting always applies to xi_2
  if (xi_1$n_sample == 1 && xi_2$n_sample > 1) {
    tmp <- xi_1
    xi_1 <- xi_2
    xi_2 <- tmp
  }

  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample

  if (!(n1 == n2 || n2 == 1)) {
    stop("Sample sizes are incompatible for broadcasting.")
  }

  # Prepare components
  f1 <- xi_1$functional_list
  f2 <- xi_2$functional_list
  Z1 <- xi_1$Z
  Z2 <- xi_2$Z

  # Replicate fd and Z if needed
  if (n2 == 1) {
    f2 <- rep_fd(f2, n1)
    Z2 <- matrix(rep(Z2, each = n1), nrow = n1, byrow = TRUE)
  }

  # Compute functional inner products
  inprod_functional <- vapply(seq_len(n1), function(i) {
    sum(vapply(seq_along(f1), function(j) {
      fda::inprod(f1[[j]][i], f2[[j]][i])
    }, numeric(1)))
  }, numeric(1))

  # Compute scalar inner products
  inprod_scalar <- rowSums(Z1 * Z2)

  # Combine results
  result <- inprod_functional + inprod_scalar
  if (n1 == 1 && n2 == 1) as.numeric(result) else result
}
