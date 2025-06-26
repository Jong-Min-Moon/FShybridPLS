# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute the First PLS Component from a Hybrid Predictor
#'
#' Computes the coefficients of the first Partial Least Squares (PLS) component based on 
#' a hybrid predictor object and a response vector, using a regularized generalized eigenvalue problem.
#'
#' @param W A \code{predictor_hybrid} object containing both functional and scalar predictors.
#' @param y A numeric response vector of length equal to the number of samples in \code{W}.
#' @param L cholesky decompisition of a regularization matrix (typically positive definite) used in the generalized eigenproblem.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{xi_hat}}{The estimated first PLS component as a \code{predictor_hybrid} object.}
#'   \item{\code{E}}{The generalized eigenvalue problem matrix \code{E = inv(L) V* t(inv(L))}.}
#'   \item{\code{V_star}}{The cross-covariance matrix between predictors and response.}
#'   \item{\code{eigen_val}}{The leading eigenvalue of \code{E}.}
#' }
#'
#' @export
get_pls_comp <- function(W, y, L){

  V_star <- get_pre_corrcov(W, y)

  invL <- Matrix::chol2inv(L)
  A <- V_star %*% t(invL)       # A = V* × t(inv(L))
  E <- invL %*% A               # E = inv(L) × A

  eigen_result <- eigen(E)
  e <- eigen_result$vectors[, 1]
  if (is.complex(e)){
    print("stop")
    return("stop")
  } 


  xi_star <- t(invL) %*% e      # xi_star solves t(L) xi = e
  xi_hat <- predictor_hybrid_from_coef(format = W, coef = xi_star)

  return(list(
    xi_hat = xi_hat,
    E = E,
    V_star = V_star,
    eigen_val = eigen_result$values[1]
  ))
}
