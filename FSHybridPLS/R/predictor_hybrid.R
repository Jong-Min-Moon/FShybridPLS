# Generated from create-rhello.Rmd: do not edit by hand

#' Create a predictor_hybrid object (S3 version)
#'
#' Constructs a hybrid predictor object that stores both scalar and functional predictors,
#' along with Jacobians and metadata.
#'
#' @param Z A numeric matrix of dimension \code{n_sample × n_scalar} representing scalar predictors.
#' @param functional_list A list of functional predictors (e.g., \code{fd} objects from the \code{fda} package).
#' @param jacobian_list A list of Jacobian matrices corresponding to the functional predictors.
#' @param n_basis_list A numeric vector indicating the number of basis functions for each functional predictor.
#' @param n_sample Number of observations (samples).
#' @param n_functional Number of functional predictors.
#' @param n_scalar Number of scalar predictors.
#'
#' @return A list of class \code{predictor_hybrid}.
#' @export
predictor_hybrid <- function(Z, functional_list, jacobian_list,
                             n_basis_list, n_sample, n_functional, n_scalar) {
  stopifnot(is.matrix(Z))
  stopifnot(length(functional_list) == n_functional)
  stopifnot(length(jacobian_list) == n_functional)
  stopifnot(length(n_basis_list) == n_functional)
  
  structure(
    list(
      Z = Z,
      functional_list = functional_list,
      jacobian_list = jacobian_list,
      n_basis_list = n_basis_list,
      n_sample = n_sample,
      n_functional = n_functional,
      n_scalar = n_scalar
    ),
    class = "predictor_hybrid"
  )
}
