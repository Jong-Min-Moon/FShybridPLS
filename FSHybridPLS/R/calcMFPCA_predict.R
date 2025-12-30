# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Prediction of MFPCA Scores
#'
#' @description
#' Calculates Multivariate Functional PCA scores for new (prediction) data based on
#' the eigenvectors and normalization factors derived from the training data.
#'
#' @param N Number of observations in the prediction set.
#' @param p Number of functional elements.
#' @param M Number of multivariate components.
#' @param weights Weights for each functional element.
#' @param npc Vector containing the number of univariate basis functions.
#' @param uniBasis List containing univariate basis decompositions for the prediction data.
#' @param normFactors Normalization factors derived from training.
#' @param vectors_train Eigenvectors derived from training.
#' @param values_train Eigenvalues derived from training.
#'
#' @return A matrix of MFPCA scores for the prediction data.
#' @export
calcMFPCA_predict <- function(N, p, M, weights, npc, uniBasis,
                              normFactors_train, vectors_train, values_train)
{
  # combine all scores
  allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores.pred}

  # block vector of weights
  allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}

  Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)


  ### Calculate scores
  # Project the prediction Z matrix onto the training eigenvectors
  scores <- Z %*% vectors_train * sqrt(N-1) # see defintion of Z above!
  scores <- as.matrix(scores %*% diag(sqrt(values_train) * normFactors_train, nrow = M, ncol = M)) # normalization

  return(scores)
}

