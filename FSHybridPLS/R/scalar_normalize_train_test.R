# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Normalize Scalar Predictors (Train/Test Split)
#'
#' Standardizes scalar predictors (Z) to have mean zero and unit variance.
#' Statistics are computed from the training set.
#'
#' @param predictor_train A `predictor_hybrid` object (train).
#' @param predictor_test A `predictor_hybrid` object (test).
#'
#' @return A list containing normalized objects and statistics.
#' @export
scalar_normalize_train_test <- function(predictor_train, predictor_test) {
  stopifnot(
    inherits(predictor_train, "predictor_hybrid"),
    inherits(predictor_test, "predictor_hybrid"),
    predictor_train$n_scalar == predictor_test$n_scalar
  )

  predictor_train_normalized <- predictor_train
  predictor_test_normalized <- predictor_test

  # 1. Compute Statistics from TRAINING data
  mean_train <- colMeans(predictor_train$Z)
  sd_train <- apply(predictor_train$Z, 2, sd)
  
  # Handle zero variance columns
  sd_train[sd_train == 0] <- 1

  # 2. Normalize TRAINING set
  predictor_train_normalized$Z <- scalar_normalize(
    predictor_train$Z,
    mean_train,
    sd_train
  )

  # 3. Normalize TEST set using TRAIN stats
  predictor_test_normalized$Z <- scalar_normalize(
    predictor_test$Z,
    mean_train,
    sd_train
  )

  return(list(
    predictor_train = predictor_train_normalized,
    predictor_test = predictor_test_normalized,
    mean_train = mean_train,
    sd_train = sd_train
  ))
}
