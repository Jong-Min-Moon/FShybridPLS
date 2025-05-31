# Generated from create-rhello.Rmd: do not edit by hand

#' Compute the Gram matrix for a basis object
#'
#' @param basis A basis object of class `basisfd` from the `fda` package
#' @return A square matrix where entry (i,j) is the inner product of basis functions i and j
#' @export
compute_gram_matrix <- function(basis) {
  stopifnot(inherits(basis, "basisfd"))
  inprod(basis, basis)  # computes ∫ b_i(t) b_j(t) dt
}

