# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct a Single-Sample Predictor Hybrid Object from Coefficients
#'
#' Reconstructs a \code{predictor_hybrid} object representing one observation, using a numeric
#' coefficient vector. This alternative constructor maps the coefficients back into their functional
#' and scalar predictor representations based on the structure of a template \code{predictor_hybrid} object.
#'
#' @param format A \code{predictor_hybrid} object that provides the structure and basis information.
#' @param coef A numeric vector containing coefficients for both functional and scalar predictors.
#'
#' @return A \code{predictor_hybrid} object with updated \code{functional_list}, \code{Z}, and \code{n_sample = 1}.
#'
#'
#' @export
predictor_hybrid_from_coef <- function(format, coef){
  M <- dim(format$gram_list[[1]])[2] # number of basis functions
  K <- length(format$functional_list) # number of functional predictors
  for (ii in 1:K){
    format$functional_list[[ii]] <- fd(
      coef = as.matrix(coef[((ii-1)*M + 1):(ii*M)]),
      basisobj = format$functional_list[[ii]]$basis
    )
  }
  format$Z <- t(as.matrix(coef[(K * M + 1):length(coef)]))
  format$n_sample <- 1
  return(format)
}
