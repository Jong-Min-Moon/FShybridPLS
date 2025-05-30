# Generated from create-rhello.Rmd: do not edit by hand

#' Add two predictor_hybrid objects
#'
#' Performs element-wise addition of two predictor_hybrid objects.
#' Functional predictors are combined using `plus.fd()` and `times.fd()`.
#'
#' @param input A predictor_hybrid object.
#' @param other Another predictor_hybrid object.
#' @param alpha A numeric scalar multiplier for `other`.
#'
#' @return A new predictor_hybrid object representing the sum.
#' @export
add.predictor_hybrid <- function(input, other, alpha = 1) {
  if (input$n_functional != other$n_functional || input$n_scalar != other$n_scalar) {
    stop("Mismatch in number of predictors.")
  }
  
  # Optionally check that functional bases are the same
  # if (!is_same_basis(input, other)) stop("Different functional bases.")

  new_functional_list <- vector("list", input$n_functional)
  for (i in seq_len(input$n_functional)) {
    new_functional_list[[i]] <- plus.fd(
      input$functional_list[[i]],
      times.fd(alpha, other$functional_list[[i]])
    )
  }
  
  new_Z <- input$Z + alpha * other$Z
  
  predictor_hybrid(
    Z = new_Z,
    functional_list = new_functional_list,
    jacobian_list = input$jacobian_list,  # Assumes identical structure
    n_basis_list = input$n_basis_list,
    n_sample = input$n_sample,
    n_functional = input$n_functional,
    n_scalar = input$n_scalar
  )
}

