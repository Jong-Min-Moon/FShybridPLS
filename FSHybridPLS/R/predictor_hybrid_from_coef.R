# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Construct a Single-Sample Predictor Hybrid Object from Coefficients
#'
#' Reconstructs a \code{predictor_hybrid} object representing one observation, using a numeric
#' coefficient vector. This alternative constructor maps the coefficients back into their functional
#' and scalar predictor representations based on the structure of a template \code{predictor_hybrid} object.
#'
#' @param format A \code{predictor_hybrid} object that provides the structure and basis information.
#' @param coef A numeric vector containing coefficients for both functional and scalar predictors.
#'
#' @return A \code{predictor_hybrid} object with updated \code{functional_list}, \code{Z}, and \code{n_sample = 1}.
#' @export
predictor_hybrid_from_coef <- function(format, coef) {
  # Extract metadata from the template
  basis_counts <- format$n_basis_list
  
  # Check if all functional predictors have the same number of basis functions
  if (length(unique(basis_counts)) != 1) {
    stop("All functional predictors must have the same number of basis functions (M).")
  }
  
  M <- basis_counts[1]  # Number of basis functions (identical across predictors)
  K <- format$n_functional # Number of functional predictors
  
  # Reconstruct each functional predictor
  for (ii in 1:K) {
    # Calculate indices for the current functional predictor in the flat vector
    start_idx <- (ii - 1) * M + 1
    end_idx   <- ii * M
    
    # Create a new fd object using the sliced coefficients and original basis
    format$functional_list[[ii]] <- fd(
      coef = as.matrix(coef[start_idx:end_idx]),
      basisobj = format$functional_list[[ii]]$basis
    )
  }
  
  # Reconstruct the scalar predictor part (Z)
  # The remaining coefficients after K*M belong to the scalar predictors
  scalar_start_idx <- K * M + 1
  format$Z <- t(as.matrix(coef[scalar_start_idx:length(coef)]))
  
  # Update sample count to 1 as this constructor creates a single observation
  format$n_sample <- 1
  
  return(format)
}
