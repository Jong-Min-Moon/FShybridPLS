# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Create a predictor_hybrid object
#'
#' @param Z Numeric matrix (n_sample x n_scalar).
#' @param functional_list List of 'fd' objects.
#' @param eval_point Evaluation points.
#' @param penalty_order Integer. Derivative order for roughness penalty (default 2).
#'
#' @export
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
