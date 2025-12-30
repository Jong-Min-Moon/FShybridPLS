# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Get Number of Samples from an fd Object
#'
#' Extracts the number of observations (replications) stored in a functional data object.
#'
#' @param fd_obj An object of class `fd`.
#' @return Integer representing the number of samples.
#' @export
n_sample.fd <- function(fd_obj) {
  if (!inherits(fd_obj, "fd")) {
    stop("Input must be of class 'fd'.")
  }
  # The coefficient matrix has dimensions (nbasis x n_sample)
  return(ncol(fd_obj$coefs))
}
