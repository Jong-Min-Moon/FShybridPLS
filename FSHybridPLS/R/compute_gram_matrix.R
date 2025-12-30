# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute the Gram matrix for a basis object
#'
#' This function calculates the inner product of a basis with itself to form
#' the Gram matrix.
#'
#' @param basis A basis object of class `basisfd` from the `fda` package.
#' @return A square matrix where entry (i,j) is the inner product of basis functions i and j.
#' @export
compute_gram_matrix <- function(basis) {
  # Ensure the input is a valid basis object from the fda package
  stopifnot(inherits(basis, "basisfd"))
  
  # Compute the inner product of the basis with itself
  # This calculates the integral of b_i(t) * b_j(t) using the metric defined on the basis
  inprod(basis, basis) 
}
