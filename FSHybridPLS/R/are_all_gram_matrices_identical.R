# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Check if all matrices in a list are identical (with tolerance for numeric comparison)
#'
#' This function iterates through a list of matrices and checks if all subsequent
#' matrices are numerically identical to the first matrix in the list.
#' It uses `all.equal` for robust comparison of numeric matrices, which accounts
#' for small floating-point differences.
#'
#' @param gram_list A list of matrices.
#' @param tolerance A numeric tolerance for comparing floating-point numbers.
#'                  Defaults to `sqrt(.Machine$double.eps)`.
#' @return A logical value: `TRUE` if all matrices are identical, `FALSE` otherwise.
#' @export
are_all_gram_matrices_identical <- function(gram_list, tolerance = sqrt(.Machine$double.eps)) {
  # Handle edge cases: empty list or list with a single matrix
  if (length(gram_list) == 0) {
    message("The input list is empty. Returning TRUE as there are no differences.")
    return(TRUE)
  }
  if (length(gram_list) == 1) {
    message("The input list contains only one matrix. Returning TRUE.")
    return(TRUE)
  }
  
  # Get the first matrix as the reference
  first_matrix <- gram_list[[1]]
  
  # Ensure the first element is a matrix
  if (!is.matrix(first_matrix)) {
    stop("The first element of 'gram_list' is not a matrix.")
  }
  
  # Iterate from the second matrix onwards and compare with the first
  for (i in 2:length(gram_list)) {
    current_matrix <- gram_list[[i]]
    
    # Ensure current element is a matrix
    if (!is.matrix(current_matrix)) {
      stop(paste("Element", i, "of 'gram_list' is not a matrix."))
    }
    
    # Compare the current matrix with the first matrix
    # all.equal returns TRUE if identical, or a character string describing differences
    # We convert the result to a logical TRUE/FALSE
    if (!isTRUE(all.equal(first_matrix, current_matrix, tolerance = tolerance))) {
      # If not equal, return FALSE immediately
      message(paste("Matrices at index 1 and", i, "are not identical."))
      return(FALSE)
    }
  }
  
  # If the loop completes, all matrices are identical
  return(TRUE)
}
