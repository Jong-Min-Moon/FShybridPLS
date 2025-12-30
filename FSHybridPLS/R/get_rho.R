# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute PLS Scores
#'
#' Projects the hybrid predictor data W onto the weight direction defined by d_vec
#' to calculate the PLS scores.
#'
#' @param d_vec A numeric vector of concatenated coefficients (functional basis coeffs + scalar weights).
#' @param W A `predictor_hybrid` object containing the data.
#'
#' @return A numeric vector of PLS scores (length n).
#' @export
get_rho <- function(d_vec, W) {
  K <- W$n_functional
  M <- W$n_basis_list[[1]] # Assumes same number of basis functions for all K
  n_scalar <- W$n_scalar
  n <- W$n_sample
  
  # Initialize rho vector
  rho <- rep(0, n)

  # --- Functional Part ---
  # Iterate over functional predictors to compute inner products <X^(k), gamma^(k)>
  for (j in 1:K) {
    # Extract coefficients for the j-th functional weight
    start_index <- (j - 1) * M + 1
    end_index <- j * M
    gamma_j <- d_vec[start_index:end_index]
    
    # Get data coefficients (samples x basis)
    Theta_j <- t(W$functional_list[[j]]$coefs)
    
    # Gram matrix
    B <- W$gram_list[[j]]
    
    # Update scores: rho += Theta * B * gamma
    # (samples x basis) * (basis x basis) * (basis x 1) -> (samples x 1)
    rho <- rho + Theta_j %*% B %*% gamma_j
  }

  # --- Scalar Part ---
  # Extract scalar weights (the remaining coefficients in d_vec)
  scalar_indices <- (K * M + 1):length(d_vec)
  scalar_weights <- d_vec[scalar_indices]
  
  # Update scores: rho += Z * weights
  rho <- rho + W$Z %*% scalar_weights
  
  return(as.vector(rho))
}
