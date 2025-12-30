# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Calculate Univariate Basis Expansion
#'
#' @description
#' Wrapper function that calls specific basis reconstruction methods (e.g., splines,
#' FPCA) to reconstruct functions from scores.
#'
#' @param type Decomposition type (e.g., "uFPCA", "splines1D").
#' @param scores Matrix of scores.
#' @param argvals Argument values.
#' @param functions Basis functions.
#' @param params Additional parameters.
#'
#' @return A \code{funData} object.
#' @export 
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

