# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Normalize Functional Predictors (Train/Test Split)
#'
#' Standardizes functional predictors to have mean zero and unit integrated variance.
#' Statistics (mean and SD) are computed from the training set and applied to the test set.
#'
#' @param train A `predictor_hybrid` object representing the training set.
#' @param test A `predictor_hybrid` object representing the testing set.
#'
#' @return A list containing:
#'   \item{predictor_train}{The normalized training object.}
#'   \item{predictor_test}{The normalized testing object.}
#'   \item{mean_train}{List of mean functions used for centering.}
#'   \item{deno_train}{Vector of scaling factors used.}
#' @export
curve_normalize_train_test <- function(train, test) {
  stopifnot(
    inherits(train, "predictor_hybrid"),
    inherits(test, "predictor_hybrid"),
    train$n_functional == test$n_functional
  )

  train_normalized <- train
  test_normalized <- test
  deno_train <- numeric(train$n_functional)
  mean_train <- list()

  # Iterate through each functional predictor
  for (k in seq_len(train$n_functional)) {
    # 1. Compute Statistics from TRAINING data
    train_fd <- train$functional_list[[k]]
    
    mean_train[[k]] <- fda::mean.fd(train_fd)
    sd_train <- fda::sd.fd(train_fd)
    
    # Scaling factor is the L2 norm of the standard deviation function
    deno_train[k] <- sqrt(fda::inprod(sd_train, sd_train))

    # Safety check for zero variance
    if (deno_train[k] <= 1e-10) {
      warning(paste("Functional predictor", k, "has near-zero variance. Skipping scaling."))
      deno_train[k] <- 1 
    }

    # 2. Normalize TRAINING set
    train_normalized$functional_list[[k]] <- curve_normalize(
      train_fd, 
      mean_train[[k]], 
      deno_train[k]
    )

    # 3. Normalize TEST set using TRAIN stats
    test_normalized$functional_list[[k]] <- curve_normalize(
      test$functional_list[[k]], 
      mean_train[[k]], 
      deno_train[k]
    )
  }

  return(list(
    predictor_train = train_normalized,
    predictor_test = test_normalized,
    mean_train = mean_train,
    deno_train = deno_train
  ))
}
