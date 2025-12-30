# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' PACE Wrapper for funData Objects
#'
#' @description
#' A wrapper function that applies the PACE algorithm (via \code{.PACE}) to
#' \code{funData} or \code{irregFunData} objects.
#'
#' @param funDataObject Training functional data object.
#' @param predData Prediction functional data object (optional).
#' @param nbasis Number of basis functions.
#' @param pve Proportion of variance explained.
#' @param npc Number of principal components.
#' @param makePD Logical. Enforce positive definiteness.
#' @param cov.weight.type Weighting type.
#'
#' @return A list containing the FPCA results formatted as \code{funData} objects.
#' @export 
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

