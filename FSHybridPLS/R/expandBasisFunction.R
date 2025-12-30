# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Reconstruct Functions from Scores and Basis
#'
#' @description
#' Calculates the linear combination of basis functions based on provided scores.
#' Used to reconstruct the approximated functions.
#'
#' @param scores Matrix of scores (coefficients).
#' @param argvals Argument values (domain).
#' @param functions Basis functions.
#'
#' @return A \code{funData} object containing the reconstructed functions.
#' @export 
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

