# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Principal Analysis by Conditional Expectation (Internal)
#'
#' @description
#' Implements the PACE method (Yao et al., 2005) for univariate Functional Principal
#' Component Analysis. This handles sparse or irregular functional data by smoothing
#' covariance surfaces.
#'
#' @param X Argument values (domain).
#' @param Y Matrix of functional data (observations).
#' @param Y.pred Matrix of functional data for prediction (optional).
#' @param nbasis Number of basis functions for smoothing.
#' @param pve Proportion of variance explained.
#' @param npc Number of principal components.
#' @param makePD Logical. Enforce positive definiteness.
#' @param cov.weight.type Weighting type for covariance smoothing.
#'
#' @return A list containing fitted components (mu, values, efunctions, scores, etc.).
#' @export 
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

