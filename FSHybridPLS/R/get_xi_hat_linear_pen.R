# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute Penalized PLS Weight Vector (Linear)
#'
#' Solves for the PLS weight direction xi that maximizes the covariance with y,
#' subject to roughness penalties on the functional components.
#'
#' @param W A `predictor_hybrid` object.
#' @param y A numeric vector of response values.
#' @param lambda A numeric vector of smoothing parameters.
#'
#' @return A single-sample `predictor_hybrid` object representing the weight direction xi.
#' @export
get_xi_hat_linear_pen <- function(W, y, lambda) {
  
  n <- W$n_sample
  K <- W$n_functional
  u <- gamma <- list()

  # --- Scalar Component ---
  # Compute scalar weights: v = Z^T * y
  v <- t(W$Z) %*% y
  
  # Initialize squared norm q with scalar part contribution
  # q = ||xi||^2_pen = sum(v^2) + functional_parts
  q <- sum(v^2)

  # --- Functional Components ---
  for (j in 1:K) {
    # Theta_t: coefficient matrix of the functional data (basis x samples)
    Theta_t <- W$functional_list[[j]]$coefs
    
    # B: Gram matrix for the j-th basis
    B <- W$gram_list[[j]]
    
    # Compute RHS of the linear system: u = B * Theta * y
    # This represents the unpenalized gradient term
    u[[j]] <- B %*% Theta_t %*% y

    # Linear System Matrix: R = Gram + lambda * Penalty
    R <- W$gram_list[[j]] + lambda[j] * W$gram_deriv2_list[[j]]
    
    # Solve for gamma_j: (Gram + Penalty) * gamma = Gram * Theta * y
    gamma[[j]] <- solve(R, u[[j]])

    # Update squared norm q:
    # Add contribution gamma^T * R * gamma (penalized norm of functional weight)
    q <- q + sum(gamma[[j]] * (R %*% gamma[[j]]))
  }
  
  # --- Combine and Normalize ---
  # Flatten gamma list and append scalar weights v
  d_vec <- c(do.call(c, gamma), v)
  
  # Normalize the vector to have unit penalized norm
  d_vec <- d_vec / sqrt(q)
  
  # Reconstruct the result as a predictor_hybrid object
  xi_hat <- predictor_hybrid_from_coef(format = W, coef = d_vec)
  
  return(xi_hat)
}
