# Generated from create-rhello.Rmd: do not edit by hand

#' @export
is_same_basis <- function(input, other) {
  all(
    length(input$functional_list) == length(other$functional_list),
    all(
      mapply(function(fd1, fd2) fda::is.eqbasis(fd1$basis, fd2$basis),
             input$functional_list, other$functional_list)
    )
  )
}
