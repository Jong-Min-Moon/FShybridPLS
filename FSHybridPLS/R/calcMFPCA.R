# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Internal Calculation of MFPCA
#'
#' @description
#' Implements the core MFPCA algorithm. It combines the univariate scores from
#' each functional element, applies weights, and performs PCA on the weighted
#' concatenated scores to extract the multivariate functional principal components.
#'
#' @param N Number of observations.
#' @param p Number of functional elements.
#' @param Bchol Cholesky decomposition of the basis integral matrix (or NULL if orthogonal).
#' @param M Number of multivariate components to calculate.
#' @param type Vector of character strings specifying the univariate expansion types.
#' @param weights Vector of weights for each functional element.
#' @param npc Vector containing the number of univariate basis functions for each element.
#' @param argvals List of argument values.
#' @param uniBasis List containing the univariate basis decompositions.
#' @param fit Logical. If TRUE, computes the fitted values (reconstruction).
#' @param approx.eigen Logical. If TRUE, uses approximate eigen decomposition.
#'
#' @return A list containing eigenvalues, eigenfunctions, scores, vectors, and normalization factors.
#' @export
calcMFPCA <- function(N, p, Bchol, M, type, weights, npc, argvals, uniBasis, fit = FALSE, approx.eigen = FALSE)
{
  # combine all scores
  allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores}

  # block vector of weights
  allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}

  Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)

  # check if approximation is appropriate (cf. irlba)
  if(approx.eigen & (M > min(N, sum(npc))/2))
  {
    warning("Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
    approx.eigen = FALSE
  }

  # check if non-orthonormal basis functions used and calculate PCA on scores
  if(is.null(Bchol))
  {
    if(approx.eigen)
    {
      tmpSVD <- irlba::irlba(as.matrix(Z), nv = M)

      vectors <- tmpSVD$v
      values <- tmpSVD$d[seq_len(M)]^2
    }
    else
    {
      if(sum(npc) > 1000)
        warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")

      e <- eigen(stats::cov(allScores) * outer(allWeights, allWeights, "*"))

      values <- e$values[seq_len(M)]
      vectors <- e$vectors[,seq_len(M)]
    }
  }
  else
  {
    if(approx.eigen)
    {
      tmpSVD <- irlba::irlba(as.matrix(Matrix::tcrossprod(Z, Bchol)), nv = M)

      vectors <- Matrix::crossprod(Bchol, tmpSVD$v)
      values <- tmpSVD$d[seq_len(M)]^2
    }
    else
    {
      if(sum(npc) > 1000)
        warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")

      e <- eigen(Matrix::crossprod(Bchol) %*% (stats::cov(allScores) * outer(allWeights, allWeights, "*")))

      values <- Re(e$values[seq_len(M)])
      vectors <- Re(e$vectors[,seq_len(M)])
    }
  }

  # normalization factors
  normFactors <- 1/sqrt(diag(as.matrix(Matrix::crossprod(Z %*% vectors))))

  ### Calculate scores
  scores <- Z %*% vectors * sqrt(N-1) # see defintion of Z above!
  scores <- as.matrix(scores %*% diag(sqrt(values) * normFactors, nrow = M, ncol = M)) # normalization

  ### Calculate eigenfunctions (incl. normalization)
  npcCum <- cumsum(c(0, npc)) # indices for blocks (-1)

  tmpWeights <- as.matrix(Matrix::crossprod(Z, Z %*%vectors))
  eFunctions <- foreach::foreach(j = seq_len(p)) %do% {
    univExpansion(type = type[j],
                  scores = 1/sqrt(weights[j] * values) * normFactors * t(tmpWeights[npcCum[j]+seq_len(npc[j]), , drop = FALSE]),
                  argvals = argvals[[j]],
                  functions = uniBasis[[j]]$functions,
                  params = uniBasis[[j]]$settings)
  }

  res <- list(values = values,
              functions = multiFunData(eFunctions),
              scores = scores,
              vectors = vectors,
              values = values,
              normFactors = normFactors,
              uniBasis = uniBasis,
              uniExpansions = uniExpansions
  )

  # calculate truncated Karhunen-Loeve representation (no mean here)
  if(fit)
    res$fit <- multivExpansion(multiFuns = res$functions, scores = scores)

  return(res)
}

