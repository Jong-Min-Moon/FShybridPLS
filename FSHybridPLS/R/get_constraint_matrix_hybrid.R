# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct penalized constraint matrix for hybrid predictors
#'
#' Computes the denominator matrix used in penalized estimation for hybrid predictors,
#' combining the Gram matrix, smoothing parameter matrix, and penalty matrix.
#'
#' Specifically, this function returns the matrix:
#' \deqn{J^\ast + \Lambda \ddot{J}^\ast}
#' where \eqn{J^\ast} is the block-diagonal Gram matrix,
#' \eqn{\Lambda} is the block-diagonal smoothing parameter matrix,
#' and \eqn{\ddot{J}^\ast} is the block-diagonal penalty matrix of second derivative inner products.
#'
#' @param W A `predictor_hybrid` object.
#' @param lambda A numeric vector of smoothing parameters, one for each functional predictor.
#'
#' @return A matrix representing the penalized constraint matrix used in estimation.
#' @export
get_constraint_matrix_hybrid <- function(W, lambda) {
  J_star <- get_gram_matrix_block(W)
  J_dotdot_star <- get_penalty_hybrid(W)
  lambda_mat <- get_smoothing_param_hybrid(W, lambda)

  J_star + lambda_mat %*% J_dotdot_star
}

