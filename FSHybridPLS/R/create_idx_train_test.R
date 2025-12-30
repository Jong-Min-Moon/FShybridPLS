# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Split Sample Indices into Train and Test Sets
#'
#' Randomly partitions the indices of a dataset into training and testing sets.
#'
#' @param n_sample Integer. The total number of samples.
#' @param train_ratio Numeric. The proportion of samples for the training set (0 to 1).
#'
#' @return A list with named elements `idx_train` and `idx_test`.
#' @export
create_idx_train_test <- function(n_sample, train_ratio) {
  # Calculate number of training samples
  n_train <- floor(n_sample * train_ratio)
  
  # Sample indices without replacement
  idx_train <- sort(sample(seq_len(n_sample), n_train, replace = FALSE))
  
  # The remaining indices form the test set
  idx_test <- sort(setdiff(seq_len(n_sample), idx_train))
  
  return(list(
    idx_train = idx_train,
    idx_test = idx_test
  ))
}
