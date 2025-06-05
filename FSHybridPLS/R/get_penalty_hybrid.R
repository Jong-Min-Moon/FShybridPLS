# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct the roughness penalty matrix for a predictor_hybrid object
#'
#' Returns a block diagonal matrix representing roughness penalties for each functional predictor
#' and zeros for the scalar predictors (scalar predictors are not penalized). Each functional component uses its own basis penalty matrix, computed via second derivative penalty via `getbasispenalty()` from the `fda` package
#'
#' @param W A `predictor_hybrid` object.
#'
#' @return A block diagonal penalty matrix of dimension `(sum(n_basis_list) + n_scalar)^2`.
#' @export
get_penalty_hybrid <- function(W){
  J_dotdot_star <- fda::getbasispenalty(W$functional_list[[1]]$basis)
  for (ii in 2:(W$n_functional)){
    J_dotdot_star <- Matrix::bdiag(
      J_dotdot_star,
      fda::getbasispenalty(W$functional_list[[ii]]$basis)
    )
    }
  J_dotdot_star <- Matrix::bdiag(
    J_dotdot_star,
    matrix(0, nrow = W$n_scalar, ncol = W$n_scalar)
    )
  return(J_dotdot_star)
  }
