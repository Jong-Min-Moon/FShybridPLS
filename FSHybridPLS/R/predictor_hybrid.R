# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Create a predictor_hybrid object (S3 version, automatic basis size and Gram matrix)
#'
#' Constructs a hybrid predictor object that stores both scalar and functional predictors,
#' and automatically computes Gram matrices (inner products of basis functions).
#'
#' @param Z A numeric matrix of dimension \code{n_sample x n_scalar} representing scalar predictors.
#' @param functional_list A list of functional predictors (e.g., \code{fd} objects from the \code{fda} package).
#' @param eval_point Points at which the functions are evaluated (stored as metadata).
#'
#' @return An object of class \code{predictor_hybrid}, containing scalar and functional data with Gram matrices.
#' @export
predictor_hybrid <- function(Z, functional_list, eval_point) {
  # Validate scalar input
  stopifnot(is.matrix(Z), is.numeric(Z))
  
  # Validate functional input
  stopifnot(is.list(functional_list))
  
  n_sample <- nrow(Z)
  n_scalar <- ncol(Z)
  n_functional <- length(functional_list)
  
  # Initialize lists for storing matrices and metadata
  gram_list <- gram_deriv2_list <- vector("list", n_functional)
  n_basis_list <- numeric(n_functional)
  
  # Iterate through each functional predictor to compute necessary matrices
  for (i in seq_len(n_functional)) {
    fd_i <- functional_list[[i]]
    stopifnot(inherits(fd_i, "fd"))
    
    # Ensure the number of functional samples matches the scalar samples
    n_fd_sample <- ncol(coef(fd_i))
    if (n_fd_sample != n_sample) {
      stop(sprintf("Functional predictor %d has %d replicates, but Z has %d rows.",
                   i, n_fd_sample, n_sample))
    }
    
    # Extract basis and compute inner product matrices
    basis_i <- fd_i$basis
    gram_list[[i]] <- compute_gram_matrix(basis_i)
    gram_deriv2_list[[i]] <- fda::getbasispenalty(basis_i) # Compute roughness penalty matrix
    
    n_basis_list[i] <- basis_i$nbasis
  }
  
  # Construct the S3 object
  structure(
    list(
      Z = Z,
      functional_list = functional_list,
      eval_point = eval_point,
      gram_list = gram_list,
      gram_deriv2_list = gram_deriv2_list, # Store the penalty matrices
      n_basis_list = n_basis_list,
      n_sample = n_sample,
      n_functional = n_functional,
      n_scalar = n_scalar
    ),
    class = "predictor_hybrid"
  )
}
