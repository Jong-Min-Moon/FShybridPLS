# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' End-to-End MFPCA Regression Pipeline
#'
#' Performs sample splitting, normalization, and dimension reduction via MFPCA (functional)
#' and PCA (scalar). Selects the optimal number of components by minimizing test set RMSE.
#'
#' @param W A `predictor_hybrid` object containing all predictors.
#' @param y Vector. Response variable.
#' @param train_ratio Numeric scalar between 0 and 1 for splitting data (default 0.7).
#' @param n_pc_max Integer. Max number of FPCs/PCs to consider (must be <= min(nbasis, n_scalar)).
#'
#' @return List containing the best model fit, optimal component count, and a list of all test RMSEs.
#' @export
fit.mfpca_regression <- function(W, y, train_ratio = 0.7, n_pc_max = 10) {

  # --- 1. Input Validation and Setup ---
  if (!inherits(W, "predictor_hybrid")) {
    stop("'W' must be a predictor_hybrid object.")
  }

  eval_point <- W$eval_point
  # Check basis size from the first functional predictor (assuming same basis for all)
  n_basis <- W$n_basis_list[[1]]
  n_scalar <- W$n_scalar

  # Safety check: Ensure n_pc_max does not exceed available dimensions
  # We need at least n_pc_max basis functions and n_pc_max scalar variables
  max_comp_safe <- min(n_pc_max, n_basis, n_scalar)
  if (n_pc_max > max_comp_safe) {
    warning(paste("n_pc_max adjusted from", n_pc_max, "to", max_comp_safe,
                  "due to dimensionality constraints (nbasis or n_scalar)."))
    n_pc_max <- max_comp_safe
  }
  if (n_pc_max < 1) stop("n_pc_max must be at least 1.")

  # Define univariate expansion parameters for MFPCA
  # We construct a list of settings for each functional predictor
  uniExpansions <- lapply(seq_len(W$n_functional), function(x) {
    list(type = "uFPCA", nbasis = n_basis, npc = n_pc_max)
  })

  # --- 2. Split and Normalize Data ---
  # Uses the split_and_normalize.all function defined previously
  processed_data <- split_and_normalize.all(W, y, train_ratio)

  # Extract processed components
  W_train <- processed_data$predictor_train
  W_test  <- processed_data$predictor_test
  y_train <- processed_data$response_train
  y_test  <- processed_data$response_test

  # --- 3. Format Functional Data for MFPCA ---
  # Convert 'fd' objects to 'funData' objects required by the MFPCA package
  # We use lapply to handle any number of functional predictors
  funData_train_list <- lapply(W_train$functional_list, function(fd_obj) {
    funData::fd2funData(fd_obj, eval_point)
  })
  funData_test_list <- lapply(W_test$functional_list, function(fd_obj) {
    funData::fd2funData(fd_obj, eval_point)
  })

  # Create multiFunData objects
  mFData_train <- funData::multiFunData(funData_train_list)
  mFData_test  <- funData::multiFunData(funData_test_list)

  # --- 4. Dimension Reduction ---

  # A. Functional Part: MFPCA
  # Compute max components (n_pc_max) once, then subset later
  MFPCA_fit <- MFPCA(mFData = mFData_train,
                     mFData_predict = mFData_test,
                     M = n_pc_max,
                     uniExpansions = uniExpansions)

  # Extract scores
  scores_fun_train <- MFPCA_fit$scores
  scores_fun_test  <- MFPCA_fit$scores.pred

  # B. Scalar Part: Standard PCA
  # Using prcomp on the normalized scalar matrix Z
  PCA_fit <- prcomp(W_train$Z, center = FALSE, scale. = FALSE)

  # Extract scores (first n_pc_max components)
  # Note: predicting on test set projects Z_test onto training rotation
  scores_scalar_train <- PCA_fit$x[, 1:n_pc_max, drop = FALSE]
  scores_scalar_test  <- predict(PCA_fit, newdata = W_test$Z)[, 1:n_pc_max, drop = FALSE]

  # --- 5. Model Selection (CV Loop) ---
  linear_fit_list <- list()
  test_rmse_list  <- numeric(n_pc_max)

  for (k in 1:n_pc_max) {
    # a. Subset scores for current complexity k
    X_train <- data.frame(
      scores_fun_train[, 1:k, drop = FALSE],
      scores_scalar_train[, 1:k, drop = FALSE]
    )
    X_test <- data.frame(
      scores_fun_test[, 1:k, drop = FALSE],
      scores_scalar_test[, 1:k, drop = FALSE]
    )

    # Set consistent column names to ensure predict() works correctly
    col_names <- c(paste0("FPC", 1:k), paste0("SPC", 1:k))
    colnames(X_train) <- col_names
    colnames(X_test)  <- col_names

    # Add response for training
    train_df <- cbind(response = y_train, X_train)

    # b. Fit Linear Model
    linear_fit_list[[k]] <- lm(response ~ ., data = train_df)

    # c. Prediction & Evaluation
    y_pred <- predict(linear_fit_list[[k]], newdata = X_test)
    test_rmse_list[k] <- sqrt(mean((y_test - y_pred)^2))
  }

  # --- 6. Final Results ---
  best_idx <- which.min(test_rmse_list)

  return(list(
    n_pc_best = best_idx,
    best_test_rmse = test_rmse_list[best_idx],
    test_rmse_by_n_pc = test_rmse_list,
    final_model = linear_fit_list[[best_idx]],
    all_models = linear_fit_list,
    # Return normalization stats to allow back-transformation of predictions if needed
    normalization_details = processed_data$details
  ))
}
