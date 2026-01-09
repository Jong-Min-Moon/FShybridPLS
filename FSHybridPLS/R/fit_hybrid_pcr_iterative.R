# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Iterative Hybrid Principal Component Regression
#'
#' @description
#' Fits a Principal Component Regression (PCR) model that integrates both scalar and 
#' functional predictors. The function performs dimensionality reduction separately 
#' for each modality and iteratively fits linear models with an increasing number 
#' of principal components.
#'
#' @details
#' The workflow consists of three main stages:
#' /enumerate{
#'   /item /strong{Scalar PCA:} Performs standard PCA on the scalar matrix /code{Z}. 
#'         Test data is projected onto the training rotation matrix.
#'   /item /strong{Functional PCA:} Performs FPCA on each functional predictor using 
#'         the /code{fda} package. Crucially, test curves are centered using the 
#'         /emph{training mean function} before projection onto training harmonics 
#'         to prevent data leakage.
#'   /item /strong{Iterative Regression:} Fits a series of linear models (/code{lm}). 
#'         Model /code{l} includes the first /code{l} principal components from 
#'         both the scalar and functional sets.
#' }
#'
#' @param W_train A /code{predictor_hybrid} object containing training data.
#' @param y_train A numeric vector of training response values.
#' @param W_test A /code{predictor_hybrid} object containing test data.
#' @param y_test A numeric vector of test response values.
#' @param n_comp_max Integer. The maximum number of principal components to include 
#'        in the regression models. Defaults to 10.
#'
#' @return A list containing:
#' /item{validation_rmse}{A numeric vector of length /code{n_comp_max} containing 
#' the RMSE on the test set for each model complexity level.}
#'
#' @importFrom stats prcomp predict lm
#' @importFrom fda pca.fd fd inprod
#' @export
fit_hybrid_pcr_iterative <- function(W_train, y_train, W_test, y_test, n_comp_max = 10) {
  
  validation_rmse <- numeric(n_comp_max)
  
  # 1. Scalar PCA (Train)
  # Computes rotation matrix on training data only
  k_z_max <- min(n_comp_max, ncol(W_train$Z))
  pca_z <- prcomp(W_train$Z, center = TRUE, scale. = TRUE)
  
  # Extract training scores and project test data onto training rotation
  all_z_scores_train <- pca_z$x[, 1:k_z_max, drop = FALSE]
  all_z_scores_test  <- predict(pca_z, newdata = W_test$Z)[, 1:k_z_max, drop = FALSE]
  
  # 2. Functional PCA (Train)
  n_funcs <- length(W_train$functional_list)
  all_f_scores_train_list <- list()
  all_f_scores_test_list  <- list()
  
  for(i in 1:n_funcs) {
    # Fit FPCA on Training Data to get harmonics and mean
    pca_fd <- fda::pca.fd(W_train$functional_list[[i]], nharm = n_comp_max)
    all_f_scores_train_list[[i]] <- pca_fd$scores
    
    # Project Test Data (Robust Manual Centering)
    mean_fd <- pca_fd$meanfd
    harmonics <- pca_fd$harmonics
    
    # Extract test coefficients and training mean coefficients
    coefs_test <- W_test$functional_list[[i]]$coefs
    coefs_mean <- as.vector(mean_fd$coefs)
    
    # Sweep subtract: Test_Centered = Test_Raw - Mean_Train
    coefs_centered <- sweep(coefs_test, 1, coefs_mean, "-")
    centered_test  <- fd(coefs_centered, W_test$functional_list[[i]]$basis)
    
    # Project centered test data onto training harmonics
    all_f_scores_test_list[[i]] <- fda::inprod(centered_test, harmonics)
  }
  
  # 3. Iterative Regression Loop
  for (l in 1:n_comp_max) {
    # Scalar Scores Subset
    k_curr <- min(l, k_z_max)
    z_train_curr <- all_z_scores_train[, 1:k_curr, drop=FALSE]
    z_test_curr  <- all_z_scores_test[, 1:k_curr, drop=FALSE]
    
    # Functional Scores Subset
    f_train_curr_list <- list()
    f_test_curr_list  <- list()
    for(i in 1:n_funcs) {
      f_train_curr_list[[i]] <- all_f_scores_train_list[[i]][, 1:l, drop=FALSE]
      f_test_curr_list[[i]]  <- all_f_scores_test_list[[i]][, 1:l, drop=FALSE]
    }
    
    # Construct Design Matrices
    df_train <- data.frame(z_train_curr, do.call(cbind, f_train_curr_list))
    df_test  <- data.frame(z_test_curr,  do.call(cbind, f_test_curr_list))
    
    # Fit Linear Model
    df_train$Y <- y_train
    model <- lm(Y ~ ., data = df_train)
    
    # Predict and Compute RMSE
    y_pred <- predict(model, newdata = df_test)
    validation_rmse[l] <- sqrt(mean((y_test - y_pred)^2))
  }
  
  return(list(validation_rmse = validation_rmse))
}
