# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Calculate Functional PCA Basis for Training and Prediction Data
#'
#' @description
#' A wrapper around the PACE (Principal Analysis by Conditional Expectation) algorithm.
#' It computes FPCA scores and eigenfunctions for the training data and simultaneously
#' computes projected scores for a prediction dataset.
#'
#' @param funDataObject The training functional data object.
#' @param funDataObject_predict The prediction functional data object.
#' @param nbasis Number of basis functions to use for smoothing.
#' @param pve Proportion of Variance Explained.
#' @param npc Number of Principal Components.
#' @param makePD Logical. Enforce positive definiteness of covariance.
#' @param cov.weight.type Weighting scheme for covariance estimation.
#'
#' @return A list containing training scores, prediction scores, orthogonality status, functions, and the mean.
#' @export  
################# modified by Jongmin Mun #####################
fpcaBasis <- function(funDataObject, funDataObject_predict, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  # Calculate FPCA for training data
  FPCA <- PACE(funDataObject, predData = NULL, nbasis, pve, npc, makePD, cov.weight.type)
  # Calculate FPCA scores for prediction data using the training parameters (implicit in PACE call structure)
  FPCA_pred <- PACE(funDataObject = funDataObject, predData = funDataObject_predict, nbasis, pve, npc, makePD, cov.weight.type)

  return(list(scores = FPCA$scores,
              scores.pred = FPCA_pred$scores,
              ortho = TRUE,
              functions = FPCA$functions,
              meanFunction = FPCA$mu
  ))
}

