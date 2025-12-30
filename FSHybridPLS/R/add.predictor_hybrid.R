# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Add two predictor_hybrid objects
#'
#' Performs element-wise addition of two \code{predictor_hybrid} objects. 
#' Supports broadcasting if one object is a single sample.
#'
#' @param xi_1 A \code{predictor_hybrid} object.
#' @param xi_2 Another \code{predictor_hybrid} object to be added.
#' @param alpha A scalar multiplier applied to \code{xi_2} before addition (default is 1).
#'
#' @return A new \code{predictor_hybrid} object representing the result of the addition.
#' @export
add.predictor_hybrid <- function(xi_1, xi_2, alpha = 1) {
  # Safe access to is.eqbasis() from the fda namespace
  is_eqbasis <- getFromNamespace("is.eqbasis", "fda")
  
  # Type checks
  if (!inherits(xi_1, "predictor_hybrid") || !inherits(xi_2, "predictor_hybrid")) {
    stop("Both inputs must be of class 'predictor_hybrid'.")
  }
  
  # Structural checks: Number of functional and scalar predictors must match
  if (xi_1$n_functional != xi_2$n_functional) {
    stop("Mismatch in number of functional predictors.")
  }
  if (xi_1$n_scalar != xi_2$n_scalar) {
    stop("Mismatch in number of scalar predictors.")
  }
  
  # Basis checks: Ensure corresponding functional predictors share the same basis
  for (i in seq_len(xi_1$n_functional)) {
    if (!is_eqbasis(xi_1$functional_list[[i]]$basis, xi_2$functional_list[[i]]$basis)) {
      stop("Functional predictors must have the same basis.")
    }
  }
  
  # Broadcasting Logic:
  # Ensure xi_1 is always the larger object (or equal) to simplify broadcasting
  if (xi_1$n_sample == 1 && xi_2$n_sample > 1) {
    tmp <- xi_1
    xi_1 <- xi_2
    xi_2 <- tmp
  }
  
  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample
  
  # Validate sample compatibility
  if (!(n1 == n2 || n2 == 1)) {
    stop("Sample sizes are incompatible for broadcasting.")
  }
  
  # Prepare components
  f1 <- xi_1$functional_list
  f2 <- xi_2$functional_list
  Z1 <- xi_1$Z
  Z2 <- xi_2$Z
  
  # Replicate fd and Z if broadcasting is needed (n2 == 1)
  if (n2 == 1 && n1 > 1) {
    f2 <- rep_fd(f2, n1)
    Z2 <- matrix(rep(c(Z2), n1), nrow = n1, byrow = TRUE)
  }
  
  # Combine functional predictors: f1 + alpha * f2
  new_functional_list <- Map(
    function(fd1, fd2) plus.fd(fd1, times.fd(alpha, fd2)),
    f1,
    f2
  )
  
  # Combine scalar predictors: Z1 + alpha * Z2
  new_Z <- Z1 + alpha * Z2
  
  # Construct and return the new predictor_hybrid object
  # Note: eval_point is carried over from the first object
  predictor_hybrid(Z = new_Z, functional_list = new_functional_list, eval_point = xi_1$eval_point)
}
