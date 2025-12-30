# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Broadcast Addition for Matrices
#'
#' Adds a single observation (vector) to every row of a matrix.
#' Checks that the input is a matrix and the operand is a single observation.
#'
#' @param input A numeric matrix.
#' @param other A numeric vector or a 1-row matrix.
#' @param alpha A scalar multiplier (default 1).
#' @return A matrix with the broadcasted addition applied.
#' @export
add_broadcast <- function(input, other, alpha = 1) {
  
  # 1. Type Check: Ensure input is a matrix
  if (!is.matrix(input)) {
    stop("Argument 'input' must be a matrix.")
  }

  # 2. Structural Check: Ensure RHS is a single observation.
  # This uses length() for vectors, and dim()[1] for matrices.
  is_matrix_other <- is.matrix(other)
  if (is_matrix_other && dim(other)[1] > 1) {
    stop("Argument 'other' must be a vector or a single-row matrix.")
  }

  # 3. Conversion: Convert a single-row matrix into a vector for broadcasting.
  # If 'other' is already a vector, this does nothing.
  if (is_matrix_other) {
    other <- as.vector(other)
  }
  
  # 4. Dimension Check: Ensure vector length matches matrix columns
  if (length(other) != ncol(input)) {
    stop("Length of 'other' must match number of columns in 'input'.")
  }

  # 5. Native Arithmetic: Use simple R arithmetic for high performance.
  # R broadcasts the vector 'other' across the rows of the matrix 'input'.
  # Transposing twice (t(t(...) ...)) ensures correct column-wise recycling 
  # logic is applied to rows.
  return(t(t(input) + alpha * other))
}
