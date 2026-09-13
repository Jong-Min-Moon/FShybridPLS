#' Split Hybrid Predictor and Response Data
#'
#' Randomly partitions the hybrid predictor object and response vector into
#' training and testing sets based on a given ratio.
#'
#' @param W_hybrid A `predictor_hybrid` object.
#' @param response A numeric (or factor) vector of length `W_hybrid$n_sample`.
#' @param train_ratio Proportion of samples for training (strictly between 0 and 1).
#'
#' @return A named list with:
#' \describe{
#'   \item{predictor_train, predictor_test}{`predictor_hybrid` objects.}
#'   \item{response_train, response_test}{Response vectors for each split.}
#'   \item{train_idx, test_idx}{Integer indices used for the split.}
#' }
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 30, n_basis = 5)
#' parts <- split_all(sim$W, sim$y, train_ratio = 0.7)
#' length(parts$response_train) > 0
#' length(parts$response_test) > 0
#' inherits(parts$predictor_train, "predictor_hybrid")
#'
#' @export
split_all <- function(W_hybrid, response, train_ratio) {
  assert_predictor_hybrid(W_hybrid, "W_hybrid")
  if (!(is.numeric(response) || is.factor(response) || is.logical(response) || is.vector(response))) {
    stop("'response' must be a vector.", call. = FALSE)
  }
  if (length(response) != W_hybrid$n_sample) {
    stop("Response length must match the number of samples in W_hybrid.", call. = FALSE)
  }
  if (!is.numeric(train_ratio) || length(train_ratio) != 1 ||
      train_ratio <= 0 || train_ratio >= 1) {
    stop("'train_ratio' must be a numeric value strictly between 0 and 1.", call. = FALSE)
  }

  # Match original split.all() behavior used in the paper replication code
  N <- W_hybrid$n_sample
  N_train <- floor(N * train_ratio)
  train_idx <- sort(sample(seq_len(N), N_train))
  test_idx <- sort(setdiff(seq_len(N), train_idx))

  Z_train <- W_hybrid$Z[train_idx, , drop = FALSE]
  Z_test <- W_hybrid$Z[test_idx, , drop = FALSE]
  fun_list_train <- lapply(W_hybrid$functional_list, function(fd_obj) fd_obj[train_idx])
  fun_list_test <- lapply(W_hybrid$functional_list, function(fd_obj) fd_obj[test_idx])

  W_train <- predictor_hybrid(
    Z = Z_train,
    functional_list = fun_list_train,
    eval_point = W_hybrid$eval_point
  )
  W_test <- predictor_hybrid(
    Z = Z_test,
    functional_list = fun_list_test,
    eval_point = W_hybrid$eval_point
  )

  list(
    predictor_train = W_train,
    predictor_test = W_test,
    response_train = response[train_idx],
    response_test = response[test_idx],
    train_idx = train_idx,
    test_idx = test_idx
  )
}
