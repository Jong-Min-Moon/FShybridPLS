# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Exact Hybrid Inner Product
#' Computes <W, xi> using basis inner product matrices.
#' @param W The hybrid predictor object (contains coefs for N samples)
#' @param xi The hybrid weight object (contains coefs for 1 direction)
#' @return A vector of scores (length N)
get_rho <- function(W, xi) {
  
  # 1. Scalar Part (Standard dot product)
  # (N x p_scalar) %*% (p_scalar x 1)
  scores <- as.matrix(W$Z) %*% as.matrix(xi$Z)
  
  # 2. Functional Part (Basis Matrix Algebra)
  # Loop through each functional variable
  for(k in seq_along(W$functional_list)) {
    
    # Extract basis and coefficients
    # Coefs_W is (Basis_Dim x N)
    # Coefs_xi is (Basis_Dim x 1)
    fd_W  <- W$functional_list[[k]]
    fd_xi <- xi$functional_list[[k]]
    
    # A. Compute J Matrix (Basis Inner Product)
    # J[i,j] = Integral(phi_i(t) * phi_j(t) dt)
    # Note: computed once per variable, very fast for B-splines
    J <- fda::inprod(fd_W$basis, fd_W$basis)
    
    # B. Compute Integral via Matrix Multiplication
    # Formula: Score = Coefs_W^T * J * Coefs_xi
    # Dimensions: (N x B) * (B x B) * (B x 1) -> (N x 1)
    
    # Optimization: Pre-multiply J * xi to get a "weighted weight"
    weighted_xi <- J %*% fd_xi$coefs 
    
    # Project data onto this weighted weight
    scores <- scores + t(fd_W$coefs) %*% weighted_xi
  }
  
  return(as.vector(scores))
}
