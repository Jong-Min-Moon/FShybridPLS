# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct block-diagonal smoothing parameter matrix
#'
#' Generates a block-diagonal matrix with smoothing parameters applied to each 
#' functional component. Each block is a scaled identity matrix, where the scaling 
#' factor corresponds to the regularization parameter for that functional component. 
#' The scalar components are not penalized and thus contribute a zero matrix block.
#'
#' @param W A `predictor_hybrid` object.
#' @param lambda A numeric vector of length equal to the number of functional components (`W$n_functional`), containing the smoothing parameters for each functional predictor.
#'
#' @return A block-diagonal matrix of size `(total_dim × total_dim)`, where the top-left blocks are scaled identity matrices for functional predictors and the bottom-right block is a zero matrix for scalar covariates.
#' @export
get_smoothing_param_hybrid <- function(W, lambda) {
  if (!inherits(W, "predictor_hybrid")) {
    stop("Input W must be of class 'predictor_hybrid'.")
  }

  if (length(lambda) != W$n_functional) {
    stop("Length of lambda must match the number of functional predictors.")
  }

  lambda_blocks <- lapply(seq_len(W$n_functional), function(ii) {
    nb <- W$functional_list[[ii]]$basis$nbasis
    lambda[ii] * diag(nb)
  })

  lambda_blocks[[W$n_functional + 1]] <- matrix(0, nrow = W$n_scalar, ncol = W$n_scalar)

  Matrix::bdiag(lambda_blocks)
}
