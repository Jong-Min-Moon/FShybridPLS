# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct block-diagonal smoothing parameter matrix
#'
#' Generates a block-diagonal matrix containing regularization parameters. 
#' Functional coefficients are penalized by scaled identity matrices, while 
#' scalar coefficients receive zero penalty.
#'
#' @param W A `predictor_hybrid` object.
#' @param lambda A numeric vector of smoothing parameters, one for each functional predictor.
#'
#' @return A sparse block-diagonal matrix. The top-left blocks are `lambda[k] * I`, 
#' and the bottom-right block is a zero matrix for scalar covariates.
#' @export
get_smoothing_param_hybrid <- function(W, lambda) {
  # Input validation
  if (!inherits(W, "predictor_hybrid")) {
    stop("Input W must be of class 'predictor_hybrid'.")
  }
  if (length(lambda) != W$n_functional) {
    stop("Length of lambda must match the number of functional predictors.")
  }

  # 1. Create penalty blocks for functional predictors
  #    For each predictor k, create a scaled identity matrix: lambda[k] * I_Mk
  lambda_blocks <- lapply(seq_len(W$n_functional), function(ii) {
    nb <- W$functional_list[[ii]]$basis$nbasis
    lambda[ii] * diag(nb)
  })

  # 2. Create zero block for scalar predictors
  #    Scalar predictors are not penalized in this matrix, so we append a p x p zero matrix.
  lambda_blocks[[W$n_functional + 1]] <- matrix(0, nrow = W$n_scalar, ncol = W$n_scalar)

  # 3. Construct the sparse block-diagonal matrix
  Matrix::bdiag(lambda_blocks)
}
