# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Normalize a scalar predictor matrix
#'
#' Centers columns by 'mean' and scales them by 'deno'.
#'
#' @param scalar_predictor Numeric matrix.
#' @param mean Numeric vector of column means.
#' @param deno Numeric vector of column standard deviations.
#' @return Normalized matrix.
#' @export
scalar_normalize <- function(scalar_predictor, mean, deno) {
  # Subtract mean from each row (broadcasted)
  # Uses the helper subtr_broadcast function defined previously
  scalar_predictor_normalized <- subtr_broadcast(scalar_predictor, mean)
  
  # Create a matrix of denominators matching the dimensions of the predictor
  deno_mat <- matrix(rep(deno, each = nrow(scalar_predictor)), 
                     nrow = nrow(scalar_predictor), 
                     ncol = length(deno))
  
  # Element-wise division
  scalar_predictor_normalized <- scalar_predictor_normalized / deno_mat
  return(scalar_predictor_normalized)
}
