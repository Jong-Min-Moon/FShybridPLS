# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Replicate a list of single-sample fd objects multiple times
#'
#' This function broadcasts a list of single-sample functional predictors
#' by replicating their coefficient columns.
#'
#' @param fd_list A list of `fd` objects, each representing a single sample.
#' @param n The number of replications desired.
#'
#' @return A list of `fd` objects, each with `n` replications.
#' @export
rep_fd <- function(fd_list, n) {
  # Validate input: must be a list of 'fd' objects
  if (!is.list(fd_list) || any(!vapply(fd_list, inherits, logical(1), "fd"))) {
    stop("Input must be a list of 'fd' objects.")
  }
  
  lapply(fd_list, function(fd_obj) {
    coef_mat <- fd_obj$coefs
    
    # Ensure the object represents a single sample (one column of coefficients)
    if (is.null(dim(coef_mat)) || ncol(coef_mat) != 1) {
      stop("Each fd object in the list must have one column of coefficients.")
    }
    
    # Create a new coefficient matrix by repeating the column 'n' times
    new_coefs <- matrix(rep(coef_mat, n), nrow = nrow(coef_mat), ncol = n)
    
    # Return a new fd object with the replicated coefficients and original basis
    fd(coef = new_coefs, basisobj = fd_obj$basis)
  })
}
