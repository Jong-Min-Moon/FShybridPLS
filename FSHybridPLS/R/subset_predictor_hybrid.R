# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Extract a single observation from a predictor_hybrid object
#'
#' @param W A predictor_hybrid object with multiple samples.
#' @param i Integer index of the sample to extract.
#' @return A single-sample predictor_hybrid object.
subset_predictor_hybrid <- function(W, i) {
  new_Z <- matrix(W$Z[i, ], nrow = 1)
  new_functional_list <- lapply(W$functional_list, function(fdobj) {
    fd(coef = matrix(coef(fdobj)[, i], ncol = 1), basisobj = fdobj$basis)
  })
  predictor_hybrid(Z = new_Z, functional_list = new_functional_list)
}

