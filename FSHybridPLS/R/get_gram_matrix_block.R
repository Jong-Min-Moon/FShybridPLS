# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct block-diagonal Gram matrix for hybrid predictor
#'
#' Returns a block-diagonal matrix containing the Gram  matrices for
#' each functional component and an identity matrix for the scalar part.
#'
#' @param obj A `predictor_hybrid` object.
#'
#' @return A block-diagonal matrix of size `(total_dim × total_dim)` where
#' functional and scalar components are arranged in order.
#' @export
get_gram_matrix_block <- function(obj) {
  if (!inherits(obj, "predictor_hybrid")) {
    stop("Input must be of class 'predictor_hybrid'.")
  }

  gram_blocks <- c(obj$gram_list, list(diag(obj$n_scalar)))
  Matrix::bdiag(gram_blocks)
}

