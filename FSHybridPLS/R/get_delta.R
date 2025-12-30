# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute hybrid regression coefficient (delta)
#'
#' Calculates the weighted average of the predictor observations W_i,
#' weighted by the scores rho_i.
#'
#' @param W A `predictor_hybrid` object.
#' @param rho A numeric vector of PLS scores.
#'
#' @return A single-sample `predictor_hybrid` object representing delta.
#' @export
get_delta <- function(W, rho) {
  # Initialize delta with the first weighted observation: rho[1] * W[1]
  delta <- scalar_mul.predictor_hybrid(subset_predictor_hybrid(W, 1), rho[1])
  
  # Iteratively add the remaining weighted observations: delta += rho[i] * W[i]
  for (i in 2:length(rho)) {
    W_i <- subset_predictor_hybrid(W, i)
    delta <- add.predictor_hybrid(delta, W_i, rho[i])
  }
  
  # Normalize by the squared norm of rho
  delta <- scalar_mul.predictor_hybrid(delta, 1 / sum(rho * rho))
  
  return(delta)
}
