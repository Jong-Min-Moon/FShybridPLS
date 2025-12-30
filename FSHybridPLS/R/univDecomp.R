# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Univariate Decomposition
#'
#' @description
#' Calculates the univariate basis decomposition for each element of a multivariate
#' functional object.
#'
#' \strong{Modification:} Accepts \code{funDataObject_predict} to facilitate
#' the pass-through of test data to the specific decomposition method (e.g., uFPCA).
#'
#' @param type Character string specifying the decomposition method (e.g., "uFPCA").
#' @param funDataObject The training functional data.
#' @param funDataObject_predict The prediction functional data.
#' @param ... Additional arguments passed to the specific decomposition function.
#'
#' @return A list containing the basis decomposition results.
#' @export 
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
  # Add prediction object to parameters if uFPCA is selected
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

