# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Split Hybrid Predictor and Response Data
#'
#' Randomly partitions the hybrid predictor object and response vector into 
#' training and testing sets based on a given ratio. Handles any number of 
#' functional predictors.
#'
#' @param W_hybrid A `predictor_hybrid` object containing all functional and scalar data.
#' @param response A numeric vector. The response variable corresponding to the samples in W_hybrid.
#' @param train_ratio A numeric scalar between 0 and 1 indicating the proportion of data to use for training.
#'
#' @return A list containing split train and test sets for predictors and response.
#' @export
split.all <- function(W_hybrid, response, train_ratio) {
  # --- Input Validation ---
  if (!inherits(W_hybrid, "predictor_hybrid")) {
    stop("'W_hybrid' must be a predictor_hybrid object.")
  }
  if (!is.vector(response) && !is.factor(response)) {
    stop("'response' must be a vector.")
  }
  if (length(response) != W_hybrid$n_sample) {
    stop("Response length must match the number of samples in W_hybrid.")
  }
  if (!is.numeric(train_ratio) || train_ratio <= 0 || train_ratio >= 1) {
    stop("'train_ratio' must be a numeric value strictly between 0 and 1.")
  }

  N <- W_hybrid$n_sample
  N_train <- floor(N * train_ratio)

  # Generate random indices
  # set.seed() should be called by the user before this function if reproducibility is needed
  train_idx <- sort(sample(seq_len(N), N_train))
  test_idx  <- sort(setdiff(seq_len(N), train_idx))

  # --- Subset Scalar Predictors ---
  Z_train <- W_hybrid$Z[train_idx, , drop = FALSE]
  Z_test  <- W_hybrid$Z[test_idx,  , drop = FALSE]

  # --- Subset Functional Predictors ---
  # Use lapply to handle any number of functional predictors (K)
  # The '[' operator on an fd object subsets the samples (replicates)
  fun_list_train <- lapply(W_hybrid$functional_list, function(fd_obj) fd_obj[train_idx])
  fun_list_test  <- lapply(W_hybrid$functional_list, function(fd_obj) fd_obj[test_idx])

  # --- Reconstruct Hybrid Objects ---
  # We reuse the original evaluation points
  W_train <- predictor_hybrid(Z = Z_train, functional_list = fun_list_train, eval_point = W_hybrid$eval_point)
  W_test  <- predictor_hybrid(Z = Z_test,  functional_list = fun_list_test,  eval_point = W_hybrid$eval_point)

  return(list(
    predictor_train = W_train,
    predictor_test  = W_test,
    response_train  = response[train_idx],
    response_test   = response[test_idx],
    train_idx       = train_idx, # Useful to return indices for reference
    test_idx        = test_idx
  ))
}
