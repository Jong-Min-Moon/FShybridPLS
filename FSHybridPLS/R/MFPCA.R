# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Multivariate Functional Principal Component Analysis (Modified)
#'
#' @description
#' Computes Multivariate Functional Principal Component Analysis (MFPCA) for functions
#' defined on different domains.
#'
#' \strong{Modification:} This version has been modified to accept an additional input,
#' \code{mFData_predict}. This allows the function to calculate PCA scores for both
#' the training data and a test/validation set simultaneously. This bypasses the
#' limitation of the standard \code{predict.MFPCAfit} function, which often does not
#' support out-of-sample prediction for this specific implementation.
#'
#' The function utilizes univariate basis expansions for each functional element
#' and combines them via PCA on the weighted concatenated scores.
#'
#' @param mFData A \code{multiFunData} object containing the training functional data.
#' @param mFData_predict A \code{multiFunData} object containing the test/validation functional data.
#' @param M Integer. The number of multivariate functional principal components to calculate.
#' @param uniExpansions A list of the same length as \code{mFData}, specifying the univariate basis expansion method for each element (e.g., "uFPCA").
#' @param weights Numeric vector. Weights for each functional element (defaults to 1).
#' @param fit Logical. If \code{TRUE}, returns the truncated Karhunen-Loève representation.
#' @param approx.eigen Logical. If \code{TRUE}, uses \code{irlba} for approximate singular value decomposition (faster for large M).
#' @param bootstrap Logical. If \code{TRUE}, computes bootstrap confidence bands.
#' @param nBootstrap Integer. Number of bootstrap iterations (required if \code{bootstrap = TRUE}).
#' @param bootstrapAlpha Numeric. Significance level for bootstrap confidence bands.
#' @param bootstrapStrat Factor. Stratification factor for bootstrapping.
#' @param verbose Logical. If \code{TRUE}, prints progress messages.
#'
#' @return An object of class \code{MFPCAfit} containing:
#' \item{values}{Eigenvalues of the principal components.}
#' \item{functions}{The estimated multivariate functional principal component functions.}
#' \item{scores}{A matrix of scores for the \code{mFData} (training) objects.}
#' \item{scores.pred}{A matrix of scores for the \code{mFData_predict} (test) objects.}
#' \item{vectors}{The eigenvectors associated with the combined scores.}
#' \item{normFactors}{Normalization factors used during calculation.}
#' \item{meanFunction}{The estimated mean function.}
#'
#' @export 
MFPCA <- function(mFData, mFData_predict, M, uniExpansions, weights = rep(1, length(mFData)), fit = FALSE, approx.eigen = FALSE,
                  bootstrap = FALSE, nBootstrap = NULL, bootstrapAlpha = 0.05, bootstrapStrat = NULL,
                  verbose = options()$verbose)
{
  if(! inherits(mFData, "multiFunData"))
    stop("Parameter 'mFData' must be passed as a multiFunData object.")

  # number of components
  p <- length(mFData)
  # number of observations
  N <- nObs(mFData)

  if(!all(is.numeric(M), length(M) == 1, M > 0))
    stop("Parameter 'M' must be passed as a number > 0.")

  if(!(is.list(uniExpansions) & length(uniExpansions) == p))
    stop("Parameter 'uniExpansions' must be passed as a list with the same length as 'mFData'.")

  if(!(is.numeric(weights) & length(weights) == p))
    stop("Parameter 'weights' must be passed as a vector with the same length as 'mFData'.")

  if(!is.logical(fit))
    stop("Parameter 'fit' must be passed as a logical.")

  if(!is.logical(approx.eigen))
    stop("Parameter 'approx.eigen' must be passed as a logical.")

  if(!is.logical(bootstrap))
    stop("Parameter 'bootstrap' must be passed as a logical.")

  if(bootstrap)
  {
    if(is.null(nBootstrap))
      stop("Specify number of bootstrap iterations.")

    if(any(!(0 < bootstrapAlpha & bootstrapAlpha < 1)))
      stop("Significance level for bootstrap confidence bands must be in (0,1).")

    if(!is.null(bootstrapStrat))
    {
      if(!is.factor(bootstrapStrat))
        stop("bootstrapStrat must be either NULL or a factor.")

      if(length(bootstrapStrat) != nObs(mFData))
        stop("bootstrapStrat must have the same length as the number of observations in the mFData object.")
    }
  }

  if(!is.logical(verbose))
    stop("Parameter 'verbose' must be passed as a logical.")

  # dimension for each component
  dimSupp <- dimSupp(mFData)

  # get type of univariate expansions
  type <- vapply(uniExpansions, function(l){l$type}, FUN.VALUE = "")

  # de-mean functions -> coefficients are also de-meaned!
  # do not de-mean in uFPCA, as PACE gives a smooth estimate of the mean (see below)
  m <- meanFunction(mFData, na.rm = TRUE) # ignore NAs in data
  for(j in seq_len(p))
  {
    if(type[j] != "uFPCA")
      mFData[[j]] <- mFData[[j]] - m[[j]]
  }

  if(verbose)
    cat("Calculating univariate basis expansions (", format(Sys.time(), "%T"), ")\n", sep = "")



  ########### Below is modified by Jongmin Mun #############################
  # calculate univariate basis expansion for all components
  # Modification: Passes both training (data) and test (data_predict) to univDecomp
  uniBasis <- mapply(
    function(expansion, data, data_predict){
      do.call(univDecomp, c(
        list(funDataObject = data, funDataObject_predict = data_predict),
        expansion)
      )
    },
    expansion = uniExpansions,
    data = mFData,
    data_predict = mFData_predict,
    SIMPLIFY = FALSE
  )
  ########### Above is modified by Jongmin Mun #############################





  # for uFPCA: replace estimated mean in m
  for(j in seq_len(p))
  {
    if(type[j] == "uFPCA")
      m[[j]] <- uniBasis[[j]]$meanFunction
  }

  # Multivariate FPCA
  npc <- vapply(uniBasis, function(x){dim(x$scores)[2]}, FUN.VALUE = 0) # get number of univariate basis functions

  if(M > sum(npc))
  {
    M <- sum(npc)
    warning("Function MFPCA: total number of univariate basis functions is smaller than given M. M was set to ", sum(npc), ".")
  }

  # check if non-orthonormal basis functions used
  if(all(foreach::foreach(j = seq_len(p), .combine = "c")%do%{uniBasis[[j]]$ortho}))
    Bchol = NULL
  else
  {
    # Cholesky decomposition of B = block diagonal of Cholesky decompositions
    Bchol <- Matrix::bdiag(lapply(uniBasis, function(l){
      if(l$ortho)
        res <- Matrix::Diagonal(n = ncol(l$scores))
      else
        res <- Matrix::chol(l$B)

      return(res)}))
  }

  if(verbose)
    cat("Calculating MFPCA (", format(Sys.time(), "%T"), ")\n", sep = "")

  mArgvals <- if (utils::packageVersion("funData") <= "1.2") {
    getArgvals(mFData)
  } else {
    funData::argvals(mFData)
  }

  res <- calcMFPCA(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                   npc = npc, argvals = mArgvals, uniBasis = uniBasis, fit = fit, approx.eigen = approx.eigen)

  res$meanFunction <- m # return mean function, too

  names(res$functions) <- names(mFData)

  #############modified by jongmin mun


  #####################################
  if(fit)
  {
    res$fit <- m + res$fit # add mean function to fits
    names(res$fit) <- names(mFData)
  }

  # give correct names
  namesList <- lapply(mFData, names)
  if(!all(vapply(namesList, FUN = is.null, FUN.VALUE = TRUE)))
  {
    if(length(unique(namesList)) != 1)
      warning("Elements have different curve names. Use names of the first element for the results.")

    row.names(res$scores) <- namesList[[1]]

    if(fit)
      for(i in seq_len(p))
      names(res$fit[[i]]) <- namesList[[1]]
  }


  if(type[1] == "uFPCA"){
    scores.pred<-calcMFPCA_predict(N = nObs(mFData_predict),
                                   p = p,
                                   M = M,
                                   weights = weights,
                                   npc = npc,
                                   uniBasis = uniBasis,
                                   normFactors = res$normFactors,
                                   vectors_train=res$vectors,
                                   values_train=res$values
    )
  }
  res$scores.pred <- scores.pred

  class(res) <- "MFPCAfit"
  return(res)
}

