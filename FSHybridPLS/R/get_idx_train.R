# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Generate Train Indices
#'
#' Randomly generates a sorted vector of indices for a training set.
#'
#' @param n_sample Integer. Number of samples.
#' @param train_ratio Numeric. Ratio of train samples.
#' @return Vector of sorted train indices.
#' @export
get_idx_train <- function(n_sample, train_ratio) {
  n_train <- floor(n_sample * train_ratio)
  idx_train <- sort(sample(seq_len(n_sample), n_train, replace = FALSE))
  return(idx_train)
}
