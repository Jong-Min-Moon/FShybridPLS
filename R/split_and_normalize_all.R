#' Split and Preprocess Hybrid Data
#'
#' Full preprocessing pipeline: train/test split, within-modality standardization
#' for functional and scalar predictors, between-modality variance balancing,
#' and response standardization (train statistics applied to test).
#'
#' @param W_hybrid A `predictor_hybrid` object.
#' @param response Numeric response vector (length `W_hybrid$n_sample`).
#' @param train_ratio Proportion of samples used for training (strictly between 0 and 1).
#'
#' @return A list with:
#' \describe{
#'   \item{predictor_train, predictor_test}{Normalized `predictor_hybrid` objects.}
#'   \item{response_train, response_test}{Standardized responses.}
#'   \item{details}{Intermediate normalization objects and split indices.}
#' }
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 40, n_basis = 5)
#' prep <- split_and_normalize_all(sim$W, sim$y, train_ratio = 0.7)
#' fit <- fit_hybridPLS(
#'   prep$predictor_train,
#'   prep$response_train,
#'   n_iter = 2,
#'   lambda = 1e-3,
#'   validation_data = list(
#'     W_test = prep$predictor_test,
#'     y_test = prep$response_test
#'   )
#' )
#' fit$validation_rmse
#'
#' @export
#' @importFrom stats sd
split_and_normalize_all <- function(W_hybrid, response, train_ratio) {
  assert_predictor_hybrid(W_hybrid, "W_hybrid")
  if (!is.numeric(response) || length(response) != W_hybrid$n_sample) {
    stop("'response' must be numeric of length W_hybrid$n_sample.", call. = FALSE)
  }

  split_result <- split_all(W_hybrid, response, train_ratio)

  curve_norm <- curve_normalize_train_test(
    split_result$predictor_train,
    split_result$predictor_test
  )
  scalar_norm <- scalar_normalize_train_test(
    curve_norm$predictor_train,
    curve_norm$predictor_test
  )
  btwn_norm <- btwn_normalize_train_test(
    scalar_norm$predictor_train,
    scalar_norm$predictor_test
  )

  response_mean_train <- mean(split_result$response_train)
  response_sd_train <- stats::sd(split_result$response_train)
  if (response_sd_train == 0) response_sd_train <- 1

  response_train_std <- (split_result$response_train - response_mean_train) / response_sd_train
  response_test_std <- (split_result$response_test - response_mean_train) / response_sd_train

  list(
    predictor_train = btwn_norm$predictor_train,
    predictor_test = btwn_norm$predictor_test,
    response_train = response_train_std,
    response_test = response_test_std,
    details = list(
      curve_normalize_result = curve_norm,
      scalar_normalize_result = scalar_norm,
      btwn_normalize_result = btwn_norm,
      response_normalize_result = list(
        mean_train = response_mean_train,
        sd_train = response_sd_train
      ),
      split_indices = list(
        train = split_result$train_idx,
        test = split_result$test_idx
      )
    )
  )
}
