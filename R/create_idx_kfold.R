#' Generate K-Fold Cross-Validation Indices
#'
#' @param n_sample Integer. Total number of samples.
#' @param n_fold Integer. Number of folds (>= 2).
#'
#' @return A list of length `n_fold`; each element has `idx_train` and `idx_valid`.
#' @examples
#' set.seed(1)
#' folds <- create_idx_kfold(20, 4)
#' length(folds)
#' length(folds[[1]]$idx_valid)
#'
#' @export
create_idx_kfold <- function(n_sample, n_fold) {
  if (!is.numeric(n_sample) || n_sample < 2) {
    stop("'n_sample' must be >= 2.", call. = FALSE)
  }
  if (!is.numeric(n_fold) || n_fold < 2 || n_fold > n_sample) {
    stop("'n_fold' must be between 2 and n_sample.", call. = FALSE)
  }
  n_sample <- as.integer(n_sample)
  n_fold <- as.integer(n_fold)

  fold_id <- sample(rep(seq_len(n_fold), length.out = n_sample))
  lapply(seq_len(n_fold), function(f) {
    list(
      idx_train = which(fold_id != f),
      idx_valid = which(fold_id == f)
    )
  })
}
