






























# Baseline method: FPCA regression


## Modified MFPCA Function for Simultaneous Training and Prediction

The MFPCA (Multivariate Functional Principal Component Analysis) function computes multivariate functional PCA for functions defined across different (dimensional) domains.
We have modified the original implementation to admit a new input, **`mFData_predict`**, allowing the function to calculate PCA scores for **both the training data and new test (or validation) data simultaneously** in a single run. This change was necessary because the standard `predict.MFPCAfit` function does not permit computing scores for new, out-of-sample functions.
The multivariate functional principal component analysis inherently relies on a **univariate basis expansion** for each functional element ($X^{(j)}$). To simplify the code modification, we use only the univariate functional PCA option and have removed other options provided by the original MFPCA function. We modifies the main function and two internal functions.


### MFPCA

```{r}
#'
#' @export MFPCA
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

```


### fpcaBasis

```{r}
################# modified by Jongmin Mun #####################
fpcaBasis <- function(funDataObject, funDataObject_predict, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  FPCA <- PACE(funDataObject, predData = NULL, nbasis, pve, npc, makePD, cov.weight.type)
  FPCA_pred <- PACE(funDataObject = funDataObject, predData = funDataObject_predict, nbasis, pve, npc, makePD, cov.weight.type)
  return(list(scores = FPCA$scores,
              scores.pred = FPCA_pred$scores,
              ortho = TRUE,
              functions = FPCA$functions,
              meanFunction = FPCA$mu
  ))
}
```


### univDecomp
```{r}
univDecomp <- function(type, funDataObject, funDataObject_predict, ...)
{
  # Parameter checking
  if(is.null(type))
    stop("Parameter 'type' is missing.")
  else
  {
    if(!is.character(type))
      stop("Parameter 'type' must be a character string. See ?univDecomp for details.")
  }

  if(class(funDataObject) != "funData")
    stop("Parameter 'funDataObject' must be a funData object.")


  # get all arguments (except for function call and type)
  params <- list(...)

  # check if type and data are of correct type
  if(is.null(type))
    stop("univDecomp: must specify 'type'.")

  if(!inherits(type, "character"))
    stop("univDecomp: 'type' must be of class character.")

  if(is.null(funDataObject))
    stop("univDecomp: must specify 'funDataObject'.")

  if(class(funDataObject) != "funData")
    stop("univDecomp: 'funDataObject' must be of class funData.")

  params$funDataObject <- funDataObject # add funDataObject (-> make sure is evaluated in correct env.)

  ############### modified by Jongmin Mun ######################
  if (type == "uFPCA"){params$funDataObject_predict <- funDataObject_predict}
  ############### modified by Jongmin Mun ######################

  res <- switch(type,
                "given" = do.call(givenBasis, params),
                "uFPCA" = do.call(fpcaBasis, params),
                "UMPCA" = do.call(umpcaBasis, params),
                "FCP_TPA" = do.call(fcptpaBasis, params),
                "splines1D" = do.call(splineBasis1D, params),
                "splines1Dpen" = do.call(splineBasis1Dpen, params),
                "splines2D" = do.call(splineBasis2D, params),
                "splines2Dpen" = do.call(splineBasis2Dpen, params),
                "fda" = do.call(fdaBasis, params),
                "DCT2D" = do.call(dctBasis2D, params),
                "DCT3D" = do.call(dctBasis3D, params),
                stop("Univariate Decomposition for 'type' = ", type, " not defined!")
  )

  if(res$ortho == FALSE & is.null(res$B))
    stop("UnivDecomp: must provide integral matrix B for non-orthonormal basis functions.")

  return(res)
}
```



## Unmodified functions

### calcBasisIntegrals
```{r}

# define global variable j, used by the foreach package and confusing R CMD CHECK
globalVariables('j')

#' Utility function that calculates matrix of basis-scalar products (one dimension)
#'
#' If the element \eqn{X^{(j)}}{X^{(j)}} is expanded in basis functions \eqn{b_i^{(j)}(t),~ i = 1, \ldots, K_j}{b_i(t)},
#' this function calculates the \eqn{K_j \times K_j}{K_j  \times K_j} matrix \eqn{B^{(jj)}}{B^{(jj)}} with entries
#' \deqn{B^{(jj)}_{mn} = \int_{\mathcal{T_j}} b_m^{(j)}(t) b_n^{(j)}(t) \mathrm{d} t}.
#'
#' @section Warning: This function is implemented only for functions on one- or two-dimensional domains.
#'
#' @param basisFunctions Array of \code{npc} basis functions of dimensions \code{npc x M1} or \code{npc x M1 x M2}.
#' @param dimSupp dimension of the support of the basis functions (1 or 2)
#' @param argvals List of corresponding x-values.
#'
#' @return A matrix containing the scalar product of all combinations of basis functions (matrix \eqn{B^{(j)}})
#'
#' @seealso \code{\link{MFPCA}}, \code{\link[funData]{dimSupp}}
#'
#' @keywords internal
calcBasisIntegrals <- function(basisFunctions, dimSupp, argvals)
{
  npc <- dim(basisFunctions)[1]

  #  integral basis matrix
  B <- array(0, dim = c(npc, npc))


  if(dimSupp == 1) # one-dimensional domain
  {
    w <- funData::.intWeights(argvals[[1]])

    for(m in seq_len(npc))
    {
      for(n in seq_len(m))
        B[m, n] <- B[n, m] <- (basisFunctions[m, ]* basisFunctions[n, ])%*%w
    }
  }
  else # two-dimesional domain (otherwise function stops before!)
  {
    w1 <- t(funData::.intWeights(argvals[[1]]))
    w2 <- funData::.intWeights(argvals[[2]])

    for(m in seq_len(npc))
    {
      for(n in seq_len(m))
        B[m, n] <- B[n, m] <-  w1 %*%(basisFunctions[m, , ]* basisFunctions[n, ,])%*%w2
    }
  }

  return(B)
}
```




### calcMFPCA
```{r}
#' Internal function that implements the MFPCA algorithm for given univariate decompositions
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
```


### calcMFPCA

```{r}
#' Internal function that implements the MFPCA algorithm for given univariate decompositions
#' @export
#'
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
```

### calcMFPCA_predict

```{r}
#' calcMFPCA_predict
#' @export
#'
calcMFPCA_predict <- function(N, p, M, weights, npc, uniBasis,
                              normFactors_train, vectors_train, values_train)
{
  # combine all scores
  allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores.pred}

  # block vector of weights
  allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}

  Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)


  ### Calculate scores
  scores <- Z %*% vectors_train * sqrt(N-1) # see defintion of Z above!
  scores <- as.matrix(scores %*% diag(sqrt(values_train) * normFactors_train, nrow = M, ncol = M)) # normalization

  return(scores)
}

```


### .PACE

```{r}
.PACE <- function(X, Y, Y.pred = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  if (is.null(Y.pred))
    Y.pred = Y
  D = NCOL(Y)
  if(D != length(X)) # check if number of observation points in X & Y are identical
    stop("different number of (potential) observation points differs in X and Y!")
  I = NROW(Y)
  I.pred = NROW(Y.pred)
  d.vec = rep(X, each = I) # use given X-values for estimation of mu
  gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
  mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec = X))
  Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for (i in seq_len(I)) {
    obs.points = which(!is.na(Y[i, ]))
    cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
    cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.G0 = diag(G.0)
  diag(G.0) = NA
  row.vec = rep(X, each = D) # use given X-values
  col.vec = rep(X, D) # use given X-values
  cov.weights <- switch(cov.weight.type,
                        none = rep(1, D^2),
                        counts = as.vector(cov.count),
                        stop("cov.weight.type ", cov.weight.type, " unknown in smooth covariance estimation"))

  npc.0 = matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0)~te(row.vec, col.vec, k = nbasis), weights = cov.weights),
                                   newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
  npc.0 = (npc.0 + t(npc.0))/2
  # no extra-option (useSymm) as in fpca.sc-method
  if (makePD) { # see fpca.sc
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE,
                            do2eigen = TRUE, trace = options()$verbose)
      as.matrix(tmp$mat)
    }
  }

  ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Chapter 8)
  w <- funData::.intWeights(X, method = "trapezoidal")
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  V <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
  ###
  evalues = replace(evalues, which(evalues <= 0), 0)
  npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
  efunctions = matrix(Winvsqrt%*%eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
  evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values[seq_len(npc)]  # use correct matrix for eigenvalue problem
  cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  ### numerical integration for estimation of sigma2
  T.len <- X[D] - X[1] # total interval length
  T1.min <- min(which(X >= X[1] + 0.25*T.len)) # left bound of narrower interval T1
  T1.max <- max(which(X <= X[D] - 0.25*T.len)) # right bound of narrower interval T1
  DIAG = (diag.G0 - diag(cov.hat))[T1.min :T1.max] # function values
  # weights
  w <- funData::.intWeights(X[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(1/(X[T1.max]-X[T1.min]) * sum(DIAG*w, na.rm = TRUE), 0) #max(1/T.len * sum(DIAG*w), 0)
  ####
  D.inv = diag(1/evalues, nrow = npc, ncol = npc)
  Z = efunctions
  Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
  fit = matrix(0, nrow = I.pred, ncol = D)
  scores = matrix(NA, nrow = I.pred, ncol = npc)
  # no calculation of confidence bands, no variance matrix
  for (i.subj in seq_len(I.pred)) {
    obs.points = which(!is.na(Y.pred[i.subj, ]))
    if (sigma2 == 0 & length(obs.points) < npc) {
      stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
    }
    Zcur = matrix(Z[obs.points, ], nrow = length(obs.points),
                  ncol = dim(Z)[2])
    ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
    scores[i.subj, ] = ZtZ_sD.inv %*% crossprod(Zcur, Y.tilde[i.subj, obs.points])
    fit[i.subj, ] = t(as.matrix(mu)) + tcrossprod(scores[i.subj, ], efunctions)
  }
  ret.objects = c("fit", "scores", "mu", "efunctions", "evalues",
                  "npc", "sigma2") # add sigma2 to output
  ret = lapply(seq_len(length(ret.objects)), function(u) get(ret.objects[u]))
  names(ret) = ret.objects
  ret$estVar <- diag(cov.hat)
  return(ret)
}
```

### PACE

```{r}
PACE <- function(funDataObject, predData = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
{
  # check inputs
  if(! class(funDataObject) %in% c("funData", "irregFunData"))
    stop("Parameter 'funDataObject' must be a funData or irregFunData object.")
  if(dimSupp(funDataObject) != 1)
    stop("PACE: Implemented only for funData objects with one-dimensional support.")
  if(methods::is(funDataObject, "irregFunData")) # for irregular functional data, use funData representation
    funDataObject <- as.funData(funDataObject)

  if(is.null(predData))
    Y.pred = NULL # use only funDataObject
  else
  {
    if(!isTRUE(all.equal(funDataObject@argvals, predData@argvals)))
      stop("PACE: funDataObject and predData must be defined on the same domains!")

    Y.pred = predData@X
  }


  if(!all(is.numeric(nbasis), length(nbasis) == 1, nbasis > 0))
    stop("Parameter 'nbasis' must be passed as a number > 0.")

  if(!all(is.numeric(pve), length(pve) == 1, 0 <= pve, pve <= 1))
    stop("Parameter 'pve' must be passed as a number between 0 and 1.")

  if(!is.null(npc) & !all(is.numeric(npc), length(npc) == 1, npc > 0))
    stop("Parameter 'npc' must be either NULL or passed as a number > 0.")

  if(!is.logical(makePD))
    stop("Parameter 'makePD' must be passed as a logical.")

  if(!is.character(cov.weight.type))
    stop("Parameter 'cov.weight.type' must be passed as a character.")


  res <- .PACE(X = funDataObject@argvals[[1]],
               Y = funDataObject@X,
               Y.pred = Y.pred,
               nbasis = nbasis, pve = pve, npc = npc, makePD = makePD,
               cov.weight.type = cov.weight.type)

  return(list(mu = funData(funDataObject@argvals, matrix(res$mu, nrow = 1)),
              values = res$evalues,
              functions = funData(funDataObject@argvals, t(res$efunctions)),
              scores = res$scores,
              scores.pred = res$scores.pred,
              fit = funData(funDataObject@argvals, res$fit),
              npc = res$npc,
              sigma2 = res$sigma2,
              estVar = funData(funDataObject@argvals, matrix(res$estVar, nrow = 1))
  ))
}
```

### expandBasisFunction

```{r}
expandBasisFunction <- function(scores, argvals = functions@argvals, functions)
{
  if(dim(scores)[2] != nObs(functions))
    stop("expandBasisFunction: number of scores for each observation and number of eigenfunctions does not match.")

  # collapse higher-dimensional functions, multiply with scores and resize the result
  d <- dim(functions@X)
  nd <- length(d)

  if(nd == 2)
    resX <- scores %*% functions@X

  if(nd == 3)
  {
    resX <- array(NA, dim = c(dim(scores)[1], d[-1]))

    for(i in seq_len(d[2]))
      resX[,i,] <- scores %*% functions@X[,i,]
  }

  if(nd == 4)
  {
    resX <- array(NA, dim = c(dim(scores)[1], d[-1]))

    for(i in seq_len(d[2]))
      for(j in seq_len(d[3]))
        resX[,i,j,] <- scores %*% functions@X[,i,j,]
  }

  if(nd > 4) # slow solution due to aperm
  {
    resX <- aperm(plyr::aaply(.data = functions@X, .margins = 3:nd,
                              .fun = function(x,y){y %*% x}, y = scores),
                  c(nd-1,nd, seq_len((nd-2))))
    dimnames(resX) <- NULL
  }

  return( funData(argvals, resX) )
}
```


### univExpansion

```{r}
univExpansion <- function(type, scores, argvals = ifelse(!is.null(functions), functions@argvals, NULL), functions, params = NULL)
{
  # Parameter checking
  if(is.null(type))
    stop("Parameter 'type' is missing.")
  else
  {
    if(!is.character(type))
      stop("Parameter 'type' must be a character string. See ?univExpansion for details.")
  }

  if(is.null(scores))
    stop("Parameter 'scores' is missing.")
  else
  {
    if(!is.matrix(scores))
      stop("Parameter 'scores' must be passed as a matrix.")
  }

  if(is.numeric(argvals))
  {
    argvals <- list(argvals)
    warning("Parameter 'argvals' was passed as a vector and transformed to a list.")
  }

  if(is.null(functions))
  {
    if(is.null(argvals))
      stop("Must pass 'argvals' if 'functions' is NULL.")
    else
    {
      if(!is.list(argvals))
        stop("Parameter 'argvals' must be passed as a list.")
    }
  }
  else
  {
    if(class(functions) != "funData")
      stop("Parameter 'functions' must be a funData object.")

    # check interaction with other parameters
    if(nObs(functions) != NCOL(scores))
      stop("Number of scores per curve does not match the number of basis functions.")

    if(!is.null(argvals) & !isTRUE(all.equal(argvals, functions@argvals)))
      stop("The parameter 'argvals' does not match the argument values of 'functions'.")
  }

  if(!is.null(params) & !is.list(params))
    stop("The parameter 'params' must be passed as a list.")

  # start calculations
  params$scores <- scores
  params$functions <- functions

  if(is.numeric(argvals))
    argvals <- list(argvals)

  params$argvals <- argvals

  res <- switch(type,
                "given" = do.call(expandBasisFunction, params),
                "uFPCA" = do.call(expandBasisFunction, params),
                "UMPCA" = do.call(expandBasisFunction, params),
                "FCP_TPA" = do.call(expandBasisFunction, params),
                "splines1D" = do.call(splineFunction1D, params),
                "splines1Dpen" = do.call(splineFunction1D, params),
                "splines2D" = do.call(splineFunction2D, params),
                "splines2Dpen" = do.call(splineFunction2Dpen, params),
                "fda" = do.call(expandBasisFunction, params),
                "DCT2D" = do.call(dctFunction2D, params),
                "DCT3D" = do.call(dctFunction3D, params),
                "default" = do.call(expandBasisFunction, params),
                stop("Univariate Expansion for 'type' = ", type, " not defined!")
  )

  return(res)
}
```

## fit.mfpca_regression
Performs sample splitting, normalization, and fits models for every combination of functional and scalar principal components up to n_pc_max. Selects the best model based on test set RMSE.

```{r}
# --- 1. MAIN FUNCTION: fit.mfpca_regression ---

#' End-to-End MFPCA Regression Pipeline
#'
#' Performs sample splitting, normalization, cross-validation (CV) to select
#' optimal functional (FPC) and scalar (SPC) components, linear regression fitting,
#' and test set prediction.
#'
#' @param W predictor_hybrid object containing all predictors.
#' @param y Vector. Response variable.
#' @param train_ratio Numeric scalar between 0 and 1 for splitting data (default 0.7).
#' @param n_pc_max Integer. Max number of FPCs/PCs to consider (must be <= min(nbasis, n_scalar)).
#' @return List containing the best model fit, optimal component count, and a list of all test RMSEs.
#' @export
fit.mfpca_regression <- function(W, y, train_ratio = 0.7, n_pc_max = 10) {

  # --- 0. Setup and Initial Checks ---
  stopifnot(inherits(W, "predictor_hybrid"))

  # Extract setup parameters
  eval_point <- W$eval_point
  n_basis <- W$n_basis_list[[1]]
  n_scalar <- W$n_scalar

  # Adjust max components to safe bounds
  max_comp <- min(n_pc_max, n_basis, n_scalar)
  if (n_pc_max != max_comp) {
    warning(paste("n_pc_max adjusted from", n_pc_max, "to", max_comp,
                  "due to limits in nbasis or n_scalar."))
    n_pc_max <- max_comp
  }

  if (n_pc_max < 1) stop("n_pc_max must be at least 1.")

  uniExpansions <- list(
    list(type = "uFPCA", nbasis = n_basis, npc = n_pc_max),
    list(type = "uFPCA", nbasis = n_basis, npc = n_pc_max)
  )


  # --- 1. Split and Normalize Full Data ---
  split_normalize_result <- split_and_normailize.all(W, y, train_ratio)
  f1_train <- split_normalize_result$predictor_train$functional_list[[1]]
  f2_train <- split_normalize_result$predictor_train$functional_list[[2]]
  f1_test <- split_normalize_result$predictor_test$functional_list[[1]]
  f2_test <- split_normalize_result$predictor_test$functional_list[[2]]


  # Format functional data
  predictor_train <- multiFunData(fd2funData(f1_train, eval_point), fd2funData(f2_train, eval_point))
  predictor_test <- multiFunData(fd2funData(f1_test, eval_point), fd2funData(f2_test, eval_point))

  y_train_norm <- split_normalize_result$response_train
  Z_train_norm <- split_normalize_result$predictor_train$Z
  Z_test_norm <- split_normalize_result$predictor_test$Z
  y_test_norm <- split_normalize_result$response_test

  # Fit MFPCA (M=n_pc_max to get max scores)
  MFPCA_fit <- MFPCA(mFData = predictor_train, mFData_predict = predictor_test,
                     M = n_pc_max, uniExpansions = uniExpansions)

  # Fit Scalar PCA (prcomp)
  PCA_fit <- prcomp(Z_train_norm, center = FALSE, scale = FALSE)

  # Extract max scores
  FPC_score_train_max <- MFPCA_fit$scores
  FPC_score_test_max <- MFPCA_fit$scores.pred
  PC_score_train_max <- PCA_fit$x
  PC_score_test_max <- predict(PCA_fit, Z_test_norm)


  # --- 3. Iterate and Fit Models (Testing FPC=npc, PC=npc) ---

  linear_fit_list <- list()
  test_rmse_list <- numeric(n_pc_max)

  for (n_pc in 1:n_pc_max) {

    # 3a. Select required number of components
    FPC_train_select <- FPC_score_train_max[, 1:n_pc, drop = FALSE]
    FPC_test_select <- FPC_score_test_max[, 1:n_pc, drop = FALSE]
    PC_train_select <- PC_score_train_max[, 1:n_pc, drop = FALSE]
    PC_test_select <- PC_score_test_max[, 1:n_pc, drop = FALSE]

    # Create column names
    col_names_fpc <- paste0("FPC", 1:n_pc)
    col_names_pc <- paste0("PC", 1:n_pc)

    # 3b. Prepare data frames for linear model
    train_data_lm <- data.frame(
      response = y_train_norm,
      FPC_train_select,
      PC_train_select
    )
    colnames(train_data_lm)[-1] <- c(col_names_fpc, col_names_pc)

    test_data_lm <- data.frame(
      FPC_test_select,
      PC_test_select
    )
    colnames(test_data_lm) <- c(col_names_fpc, col_names_pc)

    # 3c. Fit Linear Regression
    linear_fit_list[[n_pc]] <- lm(response ~ ., data = train_data_lm)

    # 3d. Predict and Calculate RMSE on Test Set
    y_pred_test <- predict(linear_fit_list[[n_pc]], test_data_lm)
    test_rmse_list[n_pc] <- sqrt(mean((y_pred_test - y_test_norm)^2))
  }

  # --- 4. Select Best Model ---
  n_pc_best <- which.min(test_rmse_list)
  best_rmse <- test_rmse_list[n_pc_best]
  best_fit <- linear_fit_list[[n_pc_best]]

  # --- 5. Return Results ---
  return(list(
    n_pc_best = n_pc_best,
    best_test_rmse = best_rmse,
    test_rmse_by_n_pc = test_rmse_list,
    final_model = best_fit,
    model_list = linear_fit_list,
    regression_coefficients = coef(best_fit)
  ))
}
```


# data generation

## matrix normal
```{r}
trial_name <- "matrix_normal"

setwd("~")
if (
  substr(getwd(), nchar(getwd())-10+1, nchar(getwd())) == "/Documents"
){
  wd_now <-  substr(getwd(), 1, nchar(getwd())-10)
}else{
  wd_now <- getwd()
}
packagepath <- paste0( wd_now, "/Documents/GitHub/fsPLS") #mac
datapath <- paste(
  wd_now,
  paste("Documents/GitHub/fsPLS/simul/data_generation", trial_name, sep = "/"),
  sep = "/"
)

devtools::load_all(packagepath)

##### parameter settings by user ######
n_rep    <- 100
n_sample <- 100
n_eval            <- 99
n_predictor_scalar <- 16
test.ratio <- 0.3
train.ratio <- 1 - test.ratio
# parameters for predictor generation
ar_slope <- 0.9
gaussian_bandwidth <- 15
freq_1 <- 10
freq_2 <- 20
signal_noise_ratio <- 0.9
scalar_signal_ratio <- 0.1

n_basis_predictor <- 20

# parameters for coefficient generation
n_basis_coef      <- 30
cycle_1 <- 5
cycle_2 <- 1

reg_coef_1_formula <- function(x){2 * sin(0.5*cycle_1*pi*x) + 4* sin(1.5*cycle_1*pi*x) + 5 * sin( 2.5 * cycle_1*pi * x)}
reg_coef_2_formula <- function(x){
  exp(
    -(cycle_2*(x-0.2)^2)/(0.04^2) )-
    exp( -(cycle_2*(x-0.4)^2)/(0.04^2) ) +
    exp( -(cycle_2*(x-0.6)^2)/(0.04^2) ) -
    exp(  -(cycle_2*(x-0.8)^2)/(0.04^2) )}



##########################################
##########################################

# predictor generating parameters
eval_point      <- seq(0,1, length = n_eval)
basis_predictor <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_predictor)
lambda_pls_base <- max(abs(inprod(basis_predictor, basis_predictor)))/max(abs(getbasispenalty(basis_predictor)))

## functional regression coefficient
basis_reg_coef  <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_coef)
reg_coef_1 <- fda::Data2fd(argvals = eval_point, y = as.vector(reg_coef_1_formula(eval_point)), basisobj = basis_reg_coef)
reg_coef_2 <- fda::Data2fd(argvals = eval_point, y = as.vector(reg_coef_2_formula(eval_point)), basisobj = basis_reg_coef)
#
## scalar regression coefficient
set.seed(1); scalar_reg_coef <- round( rnorm(n_predictor_scalar, 0, 1), 1 )







# part II. data generation
data_fun_1 <- data_fun_2<- data_scalar <- signal_functional <- signal_scalar <- signal_total <- response <- list()
U_half <- matrix(rnorm(9), nrow=3)
#U = U_half %*% t(U_half)
#U
U <- matrix(
  c(0.3585408, 0.3821649, 0.4974333,
    0.3821649, 5.6928664, 0.2608125,
    0.4974333, 0.2608125, 1.6703254
  ), nrow = 3)
for (rep_number in 1:n_rep){
  set.seed(rep_number*10)
  vec_1 <- matrix(NA, nrow=n_sample, ncol = n_eval)
  vec_2 <- matrix(NA, nrow=n_sample, ncol = n_eval)
  vec_3 <- matrix(NA, nrow=n_sample, ncol = n_eval)
  for (obs in 1:n_sample){
    realization <- rmatrixnorm(1,
                               mean = matrix(0, nrow = 3, ncol = n_eval),
                               U = U, V=solve(ar1_cov(n_eval,ar_slope,1))
    )
    vec_1[obs,] <- realization[1,]
    vec_2[obs,] <- realization[2,]
    vec_3[obs,] <- realization[3,]
  }

  data_fun_1[[rep_number]] <- fda::Data2fd(argvals = eval_point,
                                           y = t(vec_1),
                                           basisobj = basis_predictor)
  data_fun_2[[rep_number]] <- fda::Data2fd(argvals = eval_point,
                                           y = t(vec_2),
                                           basisobj = basis_predictor)
  scalar.predictor <- matrix(NA, nrow = n_sample, ncol = n_predictor_scalar)
  for (ii in 1:n_sample){
    signal.mixture.chunk <-
      scalar.predictor[ii,] <- as.numeric(
        lapply(
          chunk2(
            vec_3[ii,],
            n_predictor_scalar
          ),
          mean)
      )
  }
  data_scalar[[rep_number]] <- scalar.predictor

  signal_functional[[rep_number]] <- inprod(reg_coef_1, data_fun_1[[rep_number]]) + inprod(reg_coef_2, data_fun_2[[rep_number]])
  signal_scalar[[rep_number]]  <- t(data_scalar[[rep_number]] %*% scalar_reg_coef)
} # less than 10 seconds on laptop


signal_multiplier <- (var(unlist(signal_functional)) / var(unlist(signal_scalar)))/9

for (rep_number in 1:n_rep){
  data_scalar[[rep_number]]   <- sqrt(signal_multiplier) * data_scalar[[rep_number]]
  signal_scalar[[rep_number]] <- sqrt(signal_multiplier) * signal_scalar[[rep_number]]
  signal_total[[rep_number]] <- signal_functional[[rep_number]] +  signal_scalar[[rep_number]]
}


# generate response
noise_var <- 1/9 * var(unlist(signal_total))
for (rep_number in 1:n_rep){
  response[[rep_number]] <- signal_total[[rep_number]]+ sqrt(noise_var)*rnorm(n_sample, 0, 1)
}

save(response, data_fun_1, data_fun_2,data_scalar, response, response, eval_point, n_sample, n_basis_predictor, n_rep, lambda_pls_base, noise_var,file = "simul_data.RData")
save.image(paste(datapath, paste0(trial_name, ".RData"), sep = "/")) # creating ".RData" in current working directory


```
