# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Replace a single observation in a predictor_hybrid object
#'
#' Replaces the i-th observation of a predictor_hybrid object with a new
#' single-sample predictor_hybrid object.
#'
#' @param W A predictor_hybrid object with one or more samples.
#' @param i An integer index specifying the observation to replace.
#' @param new_W A single-sample predictor_hybrid object to use for replacement.
#'
#' @return A predictor_hybrid object with the i-th observation replaced.
#' @export
replace_obs_hybrid <- function(W, i, new_W) {
  # Input validation
  if (!inherits(W, "predictor_hybrid") || !inherits(new_W, "predictor_hybrid")) {
    stop("Both W and new_W must be of class 'predictor_hybrid'.")
  }
  if (new_W$n_sample != 1) {
    stop("new_W must be a single-sample predictor_hybrid object.")
  }
  if (i < 1 || i > W$n_sample) {
    stop(paste("Index i must be between 1 and", W$n_sample))
  }
  if (W$n_scalar != new_W$n_scalar) {
    stop("Mismatch in number of scalar predictors.")
  }
  if (W$n_functional != new_W$n_functional) {
    stop("Mismatch in number of functional predictors.")
  }

  # Check for compatible basis objects
  is_eqbasis <- getFromNamespace("is.eqbasis", "fda")
  for (j in seq_len(W$n_functional)) {
    if (!is_eqbasis(W$functional_list[[j]]$basis, new_W$functional_list[[j]]$basis)) {
      stop("Functional predictors must have the same basis objects.")
    }
  }

  # Replace the i-th observation in the scalar matrix Z
  W$Z[i, ] <- new_W$Z[1, ]

  # Replace the i-th observation in each functional predictor
  for (j in seq_along(W$functional_list)) {
    # The coefficients are stored as a matrix, with columns corresponding to samples
    W$functional_list[[j]]$coefs[, i] <- new_W$functional_list[[j]]$coefs[, 1]
  }

  return(W)
}
