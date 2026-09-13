# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Create a predictor_hybrid object
#'
#' Constructs an S3 object of class \code{predictor_hybrid} that stores scalar
#' covariates and functional predictors in a joint representation, along with
#' precomputed Gram and roughness-penalty matrices used by Hybrid PLS.
#'
#' @param Z Numeric matrix (n_sample x n_scalar).
#' @param functional_list List of 'fd' objects.
#' @param eval_point Evaluation points.
#' @param penalty_order Integer. Derivative order for roughness penalty (default 2).
#'
#' @return An object of class \code{predictor_hybrid}, a named list with:
#' \describe{
#'   \item{\code{Z}}{Scalar predictor matrix (\code{n_sample} x \code{n_scalar}).}
#'   \item{\code{functional_list}}{List of \code{fd} objects for the functional predictors.}
#'   \item{\code{gram_list}}{List of Gram (inner-product) matrices for each functional basis.}
#'   \item{\code{gram_deriv_list}}{List of roughness-penalty matrices for each functional basis.}
#'   \item{\code{eval_point}}{Numeric vector of evaluation points for the functional domain.}
#'   \item{\code{n_basis_list}}{Number of basis functions for each functional predictor.}
#'   \item{\code{n_sample}}{Number of samples.}
#'   \item{\code{n_functional}}{Number of functional predictors.}
#'   \item{\code{n_scalar}}{Number of scalar predictors.}
#' }
#'
#' @export
#' @examples
#' library(fda)
#' # Create basis and functional data
#' fbasis <- create.bspline.basis(c(0, 1), 5)
#' coefs <- matrix(rnorm(15), 5, 3)
#' fd_obj <- fd(coefs, fbasis)
#'
#' # Scalar data
#' Z_mat <- matrix(rnorm(9), 3, 3)
#'
#' # Construct hybrid predictor
#' W <- predictor_hybrid(
#'   Z = Z_mat,
#'   functional_list = list(fd_obj),
#'   eval_point = seq(0, 1, 0.1)
#' )
#' class(W)
#' W$n_sample
#' W$n_functional
#' W$n_scalar
predictor_hybrid <- function(Z, functional_list, eval_point, penalty_order = 2) {
  stopifnot(is.matrix(Z), is.numeric(Z), is.list(functional_list))
  
  n_sample <- nrow(Z)
  n_functional <- length(functional_list)
  
  gram_list <- vector("list", n_functional)
  gram_deriv_list <- vector("list", n_functional)
  n_basis_list <- numeric(n_functional)
  
  for (i in seq_len(n_functional)) {
    fd_i <- functional_list[[i]]
    basis_i <- fd_i$basis
    
    # 1. Validation
    if (ncol(fd_i$coefs) != n_sample) {
      stop(sprintf("Functional predictor %d sample count mismatch.", i))
    }
    
    # 2. Pre-compute Gram Matrix (Inner product of basis functions)
    # This is the "J" matrix: J_ij = <phi_i, phi_j>
    gram_list[[i]] <- fda::eval.penalty(basis_i, Lfdobj = 0)
    
    # 3. Pre-compute Penalty Matrix
    # This is the "R" matrix: R_ij = <D^2 phi_i, D^2 phi_j>
    gram_deriv_list[[i]] <- fda::eval.penalty(basis_i, Lfdobj = penalty_order)
    
    n_basis_list[i] <- basis_i$nbasis
  }
  
  structure(
    list(
      Z = Z,
      functional_list = functional_list,
      gram_list = gram_list,            # J matrices
      gram_deriv_list = gram_deriv_list, # R matrices
      eval_point = eval_point,
      n_basis_list = n_basis_list,
      n_sample = n_sample,
      n_functional = n_functional,
      n_scalar = ncol(Z)
    ),
    class = "predictor_hybrid"
  )
}