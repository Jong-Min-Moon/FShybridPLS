# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute First PLS Component (Direct Hilbert Averaging)
#'
#' Computes the first PLS component as a normalized weighted average of predictor objects,
#' where weights are given by the centered response values.
#'
#' @param W A `predictor_hybrid` object containing `n_sample` observations.
#' @param y A numeric response vector of length `n_sample`.
#'
#' @return A normalized `predictor_hybrid` object representing the first PLS component.
#' @export
get_pls_comp_simple <- function(W, y) {
  n <- W$n_sample
  if (length(y) != n) stop("Length of y must match number of samples in W.")

  # Center the response (optional, but often improves interpretability)
  y_centered <- y - mean(y)

  # Initialize zero element
  xi_sum <- predictor_hybrid_from_coef(W, coef = rep(0, sum(W$n_basis_list) + W$n_scalar))
  
  # Accumulate weighted sum: sum_i y_i * W_i
  for (i in 1:n) {
    Wi <- subset_predictor_hybrid(W, i)
    xi_sum <- add.predictor_hybrid(xi_sum, Wi, alpha = y_centered[i])
  }

  # Normalize
  norm_val <- sqrt(inprod.predictor_hybrid(xi_sum))
  xi_hat <- add.predictor_hybrid(xi_sum, xi_sum, alpha = 1 / norm_val)

  return(xi_hat)
}

