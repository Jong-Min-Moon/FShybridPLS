# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Split and Preprocess Hybrid Data
#'
#' Executes the full preprocessing pipeline: splitting data, standardizing within
#' modalities (functional/scalar), balancing variance between modalities, and 
#' standardizing the response.
#'
#' @param W_hybrid The full `predictor_hybrid` object.
#' @param response The numeric response vector.
#' @param train_ratio Numeric scalar (0-1) for train/test split.
#'
#' @return A list containing processed datasets and normalization parameters.
#' @export
split_and_normalize.all <- function(W_hybrid, response, train_ratio) {
  
  # 1. Split Data
  # Note: Requires the previously defined split.all function
  split_result <- split.all(W_hybrid, response, train_ratio)

  # 2. Normalize Functional Predictors (Within)
  curve_norm <- curve_normalize_train_test(
    split_result$predictor_train, 
    split_result$predictor_test
  )
  
  # 3. Normalize Scalar Predictors (Within)
  scalar_norm <- scalar_normalize_train_test(
    curve_norm$predictor_train, 
    curve_norm$predictor_test
  )
  
  # 4. Normalize Between Modalities (Balance Variances)
  btwn_norm <- btwn_normalize_train_test(
    scalar_norm$predictor_train,
    scalar_norm$predictor_test
  )
  
  # 5. Normalize Response Variable
  response_mean_train <- mean(split_result$response_train)
  response_sd_train   <- sd(split_result$response_train)
  
  # Prevent division by zero if response is constant
  if (response_sd_train == 0) response_sd_train <- 1
  
  response_train_std <- (split_result$response_train - response_mean_train) / response_sd_train
  response_test_std  <- (split_result$response_test  - response_mean_train) / response_sd_train
  
  response_norm <- list(mean_train = response_mean_train, sd_train = response_sd_train)

  # Return comprehensive results
  return(list(
    predictor_train = btwn_norm$predictor_train,
    predictor_test  = btwn_norm$predictor_test,
    response_train  = response_train_std,
    response_test   = response_test_std,
    
    # Store intermediate results and statistics for reference/reconstruction
    details = list(
      curve_normalize_result  = curve_norm,
      scalar_normalize_result = scalar_norm,
      btwn_normalize_result   = btwn_norm,
      response_normalize_result = response_norm,
      split_indices = list(train = split_result$train_idx, test = split_result$test_idx)
    )
  ))
}
