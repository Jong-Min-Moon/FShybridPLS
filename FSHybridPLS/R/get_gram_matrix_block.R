# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct block-diagonal Gram matrix for hybrid predictor
#'
#' Creates a unified Gram matrix by placing the pre-computed functional Gram matrices
#' and a scalar identity matrix into a block-diagonal structure.
#'
#' @param obj A `predictor_hybrid` object.
#'
#' @return A sparse block-diagonal matrix of size `(total_dim x total_dim)`, where
#' `total_dim` is the sum of all functional basis sizes plus the number of scalar predictors.
#' @export
get_gram_matrix_block <- function(obj) {
  # Input validation
  if (!inherits(obj, "predictor_hybrid")) {
    stop("Input must be of class 'predictor_hybrid'.")
  }

  # 1. Retrieve the list of functional Gram matrices (pre-computed in the object)
  #    These correspond to the J^(k) blocks.
  gram_blocks <- obj$gram_list
  
  # 2. Append the Identity matrix for the scalar components
  #    This corresponds to the I_p block.
  #    We use 'diag' to create an Identity matrix of size n_scalar x n_scalar.
  gram_blocks[[length(gram_blocks) + 1]] <- diag(obj$n_scalar)
  
  # 3. Construct the sparse block-diagonal matrix
  #    Matrix::bdiag efficiently handles the block construction.
  Matrix::bdiag(gram_blocks)
}
