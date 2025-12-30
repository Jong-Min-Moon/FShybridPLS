# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Residualize the hybrid predictor
#'
#' Updates each observation in W by subtracting the projection onto rho (rho_i * delta).
#'
#' @param W The current `predictor_hybrid` object.
#' @param rho The current PLS score vector.
#' @param delta The hybrid regression coefficient (single-sample object).
#'
#' @return The residualized `predictor_hybrid` object.
#' @export
residualize_predictor <- function(W, rho, delta) {
  n <- length(rho)
  W_res <- W
  
  # Iterate through each sample to compute residuals
  for (i in 1:n) {
    # Extract the i-th observation
    W_i <- subset_predictor_hybrid(W, i)
    
    # Calculate residual: W_i_new = W_i - rho[i] * delta
    W_res_i <- subtr.predictor_hybrid(W_i, delta, rho[i])
    
    # Replace the observation in the result object
    W_res <- replace_obs_hybrid(W_res, i, W_res_i)
  }
  
  return(W_res)
}
