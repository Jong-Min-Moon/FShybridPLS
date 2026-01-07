# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Penalized Inner Product (Algebraic)
#'
#' @export
inprod_pen.predictor_hybrid <- function(xi_1, xi_2 = NULL, lambda) {
  # 1. Base Inner Product (Unpenalized)
  base_val <- inprod.predictor_hybrid(xi_1, xi_2)
  
  if (is.null(xi_2)) xi_2 <- xi_1
  
  # Determine dimensions for manual penalty calculation
  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample
  is_broadcasting <- (n2 == 1 && n1 > 1)
  # Handle the swap if inprod swapped them, or just rely on logic
  if (n1 == 1 && n2 > 1) {
     tmp <- xi_1; xi_1 <- xi_2; xi_2 <- tmp
     n1 <- xi_1$n_sample; n2 <- xi_2$n_sample
     is_broadcasting <- TRUE
  }

  penalty_val <- numeric(n1)

  # 2. Add Roughness Penalty
  for (k in seq_len(xi_1$n_functional)) {
    if (lambda[k] == 0) next # Skip if no penalty
    
    R  <- xi_1$gram_deriv_list[[k]] # (M x M) - The penalty matrix
    C1 <- xi_1$functional_list[[k]]$coefs
    C2 <- xi_2$functional_list[[k]]$coefs
    
    R_C2 <- R %*% C2
    
    if (is_broadcasting) {
      term <- as.vector(t(C1) %*% R_C2)
    } else {
      term <- colSums(C1 * R_C2)
    }
    
    penalty_val <- penalty_val + (lambda[k] * term)
  }
  
  return(base_val + penalty_val)
}
