# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Check if two functional data lists share the same basis
#'
#' @description
#' Checks if two composite functional data objects share the exact same basis functions
#' by comparing their internal list structures.
#'
#' @param input A list containing functional data objects (expected to have a `$functional_list` structure).
#' @param other Another list containing functional data objects to compare against.
#' @return `TRUE` if both lists have the same length and every corresponding basis is equivalent; otherwise `FALSE`.
#' @export
is_same_basis <- function(input, other) {
  all(
    # Check if both lists have the same number of functional elements
    length(input$functional_list) == length(other$functional_list),
    
    # Check if the basis for each corresponding element is identical
    all(
      mapply(function(fd1, fd2) {
        # Use fda package function to check basis equality
        fda::is.eqbasis(fd1$basis, fd2$basis)
      }, input$functional_list, other$functional_list)
    )
  )
}
