# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Inter-Modality Normalization (Between Functional and Scalar)
#'
#' Scales the scalar predictors (Z) so that their total variability is comparable
#' to the total variability of the functional predictors. This addresses the 
#' scale invariance issue in PLS by weighting the scalar part.
#'
#' @param train A `predictor_hybrid` object (train).
#' @param test A `predictor_hybrid` object (test).
#'
#' @return A list containing normalized objects and the scaling factor.
#' @export
btwn_normalize_train_test <- function(train, test) {
  stopifnot(inherits(train, "predictor_hybrid"), inherits(test, "predictor_hybrid"))

  train_normalized <- train
  test_normalized <- test
  Z_s <- train$Z

  # --- 1. Calculate Numerator (Total Variance of Functional Parts) ---
  omega_numerator <- 0
  for (k in seq_len(train$n_functional)) {
    # Compute inner product matrix <X_i, X_j>
    # The diagonal elements <X_i, X_i> represent the squared L2 norm ||X_i||^2
    functional_inprod_matrix <- fda::inprod(
      train$functional_list[[k]],
      train$functional_list[[k]]
    )
    omega_numerator <- omega_numerator + sum(diag(functional_inprod_matrix))
  }

  # --- 2. Calculate Denominator (Total Variance of Scalar Part) ---
  omega_denominator <- sum(Z_s^2)

  # --- 3. Compute Weight omega ---
  if (omega_denominator == 0) {
    warning("Scalar component has zero variance. Skipping inter-modality scaling.")
    omega <- 1
  } else {
    omega <- omega_numerator / omega_denominator
  }

  # --- 4. Apply Scaling (sqrt(omega)) ---
  # This effectively weights the scalar inner product by omega
  scaling_factor <- sqrt(omega)
  
  train_normalized$Z <- scaling_factor * train$Z
  test_normalized$Z  <- scaling_factor * test$Z

  return(list(
    predictor_train = train_normalized,
    predictor_test = test_normalized,
    scaling_factor = scaling_factor
  ))
}
