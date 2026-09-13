# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Normalize a functional data object
#'
#' Subtracts a mean function and divides by a scalar constant.
#' Manual broadcasting is used to ensure robust subtraction across samples.
#'
#' @keywords internal
#' @param functional An `fd` object to normalize.
#' @keywords internal
#' @param mean_functional An `fd` object representing the mean.
#' @keywords internal
#' @param deno A numeric scalar for scaling (division).
#'
#' @return A normalized `fd` object.
#' @importFrom stats coef
#' @importFrom fda fd times.fd
curve_normalize <- function(functional, mean_functional, deno) {
  # Replicate the mean coefficients to match the number of samples in 'functional'
  # This creates a matrix of identical mean columns for broadcasting
  mean_coefs_matrix <- coef(mean_functional) %*% matrix(1, ncol = ncol(coef(functional)))
  
  # Subtract mean and create new fd object
  functional_normalized <- fd(
    coef = coef(functional) - mean_coefs_matrix,
    basisobj = functional$basis
  )
  
  # Apply scaling factor
  # times.fd usually handles scalar multiplication well
  functional_normalized <- times.fd(1 / deno, functional_normalized)
  
  return(functional_normalized)
}