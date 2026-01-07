# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Algebraic Inner Product for Hybrid Objects
#' 
#' Uses pre-computed Gram matrices to calculate inner products 
#' without numerical integration.
#'
#' @export
inprod.predictor_hybrid <- function(xi_1, xi_2 = NULL) {
  if (is.null(xi_2)) xi_2 <- xi_1
  
  # --- 1. Identify Broadcasting Case ---
  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample
  
  # Case A: 1-to-1 (Self-norm or pair-wise)
  # Case B: N-to-1 (Broadcasting / Projection)
  is_broadcasting <- (n2 == 1 && n1 > 1)
  
  if (n1 != n2 && !is_broadcasting) {
    # Attempt to swap if the user passed (1, N) instead of (N, 1)
    if (n1 == 1 && n2 > 1) {
       tmp <- xi_1; xi_1 <- xi_2; xi_2 <- tmp
       n1 <- xi_1$n_sample; n2 <- xi_2$n_sample
       is_broadcasting <- TRUE
    } else {
       stop("Incompatible sample sizes.")
    }
  }

  # --- 2. Scalar Component ---
  # Z1: (N x P), Z2: (N x P) or (1 x P)
  if (is_broadcasting) {
    # Project rows of Z1 onto single vector Z2
    # Result: (N x 1) -> vector
    term_scalar <- as.vector(xi_1$Z %*% t(xi_2$Z)) 
  } else {
    # Element-wise row sums
    term_scalar <- rowSums(xi_1$Z * xi_2$Z)
  }
  
  # --- 3. Functional Component (Algebraic) ---
  term_functional <- numeric(n1)
  
  for (k in seq_len(xi_1$n_functional)) {
    J <- xi_1$gram_list[[k]]       # (M x M)
    C1 <- xi_1$functional_list[[k]]$coefs # (M x N)
    C2 <- xi_2$functional_list[[k]]$coefs # (M x N) or (M x 1)
    
    # Calculate weighted coefficients: J * C2
    J_C2 <- J %*% C2 # (M x N) or (M x 1)
    
    if (is_broadcasting) {
      # Proj = C1^T * J * c2
      # Dimension: (N x M) * (M x 1) = (N x 1)
      term_functional <- term_functional + as.vector(t(C1) %*% J_C2)
    } else {
      # Pairwise = diag(C1^T * J * C2)
      # Efficiently: colSums(C1 * (J * C2))
      # This calculates dot product of each column pair
      term_functional <- term_functional + colSums(C1 * J_C2)
    }
  }
  
  return(term_scalar + term_functional)
}
