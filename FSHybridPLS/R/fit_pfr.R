# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Fit Penalized Functional Regression (PFR)
#'
#' This function formats training and testing data, fits a PFR model using the `refund` package,
#' and evaluates predictive performance via RMSE. It explicitly handles the "unrolling" of
#' scalar matrix predictors to avoid common dimension errors in `gam`/`pfr`.
#'
#' @param W_train A list containing training predictors:
#'   \itemize{
#'     \item \code{Z}: Matrix of scalar predictors.
#'     \item \code{functional_list}: List of `fd` objects (functional predictors).
#'     \item \code{eval_point}: Vector of evaluation time points.
#'   }
#' @param y_train Vector of scalar response values for training.
#' @param W_test A list containing testing predictors (same structure as \code{W_train}).
#' @param y_test Vector of scalar response values for testing (used for RMSE calculation).
#'
#' @return A numeric value representing the Root Mean Squared Error (RMSE) on the test set.
#'         Returns \code{NA} if the model fit fails.
#' @export
fit_pfr <- function(W_train, y_train, W_test, y_test) {
  
  # Extract common evaluation points for the functional curves
  eval_pts <- W_train$eval_point
  
  # ------------------------------------------------------------------
  # 1. Prepare Scalar Predictors (Unrolling)
  # ------------------------------------------------------------------
  # Convert the scalar matrix Z into a data frame with explicit column names 
  # (Z1, Z2, ...). This "unrolling" is critical because `pfr`/`gam` can 
  # throw "variable length differs" errors if predictors are stored as a 
  # single matrix column inside the data frame during prediction.
  Z_train_df <- as.data.frame(W_train$Z)
  colnames(Z_train_df) <- paste0("Z", 1:ncol(Z_train_df))
  
  Z_test_df  <- as.data.frame(W_test$Z)
  colnames(Z_test_df) <- paste0("Z", 1:ncol(Z_test_df))
  
  # ------------------------------------------------------------------
  # 2. Construct Data Frames for Model Fitting
  # ------------------------------------------------------------------
  # Combine response, unrolled scalars, and evaluated functional matrices.
  
  # Training Data
  df_train <- cbind(response = y_train, Z_train_df)
  
  # Evaluate functional objects (fd) at discrete points to create predictor matrices.
  # t() is used because eval.fd returns (time x subjects), but pfr expects (subjects x time).
  df_train$F1 <- t(eval.fd(eval_pts, W_train$functional_list[[1]])) 
  df_train$F2 <- t(eval.fd(eval_pts, W_train$functional_list[[2]]))
  
  # Testing Data
  df_test <- cbind(response = y_test, Z_test_df)
  df_test$F1 <- t(eval.fd(eval_pts, W_test$functional_list[[1]]))
  df_test$F2 <- t(eval.fd(eval_pts, W_test$functional_list[[2]]))
  
  # ------------------------------------------------------------------
  # 3. Dynamically Build the Formula
  # ------------------------------------------------------------------
  # Create the scalar part of the formula: "Z1 + Z2 + ... + Zp"
  scalar_formula_part <- paste(colnames(Z_train_df), collapse = " + ")
  
  # Construct the full formula string.
  # lf(): Defines a linear functional term.
  # presmooth = 'fpca.sc': Applies FPCA to smooth functions and reduce dimensionality 
  # prior to regression, which is robust for noisy data.
  f_str <- paste0(
    "response ~ ", scalar_formula_part, " + ",
    "lf(F1, argvals = eval_pts, presmooth = 'fpca.sc') + ",
    "lf(F2, argvals = eval_pts, presmooth = 'fpca.sc')"
  )
  
  # ------------------------------------------------------------------
  # 4. Fit Model and Predict
  # ------------------------------------------------------------------
  tryCatch({
    # Fit the PFR model.
    # method = "GCV.Cp": Uses Generalized Cross-Validation for smoothing parameter selection.
    # gamma = 1.2: Inflation factor for GCV to prevent overfitting (smoother results).
    fit <- pfr(as.formula(f_str), data = df_train, method = "GCV.Cp", gamma = 1.2)
    
    # Predict on test data
    pred <- predict(fit, newdata = df_test)
    
    # Calculate RMSE
    rmse <- sqrt(mean((y_test - pred)^2))
    return(rmse)
    
  }, error = function(e) {
    # Graceful error handling prevents the entire simulation loop from crashing
    message("PFR Error: ", e$message)
    return(NA)
  })
}
