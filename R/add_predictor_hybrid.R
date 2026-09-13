#' Add two predictor_hybrid objects
#'
#' Performs element-wise addition of two \code{predictor_hybrid} objects.
#' Supports broadcasting if one object is a single sample.
#'
#' @param xi_1 A \code{predictor_hybrid} object.
#' @param xi_2 Another \code{predictor_hybrid} object to be added.
#' @param alpha A scalar multiplier applied to \code{xi_2} before addition (default is 1).
#'
#' @return A new \code{predictor_hybrid} object representing the result of the addition.
#' @keywords internal
#' @importFrom fda plus.fd times.fd
add_predictor_hybrid <- function(xi_1, xi_2, alpha = 1) {
  if (!inherits(xi_1, "predictor_hybrid") || !inherits(xi_2, "predictor_hybrid")) {
    stop("Both inputs must be of class 'predictor_hybrid'.", call. = FALSE)
  }
  if (xi_1$n_functional != xi_2$n_functional) {
    stop("Mismatch in number of functional predictors.", call. = FALSE)
  }
  if (xi_1$n_scalar != xi_2$n_scalar) {
    stop("Mismatch in number of scalar predictors.", call. = FALSE)
  }

  for (i in seq_len(xi_1$n_functional)) {
    if (!is_same_basis(xi_1$functional_list[[i]]$basis, xi_2$functional_list[[i]]$basis)) {
      stop("Functional predictors must have the same basis.", call. = FALSE)
    }
  }

  if (xi_1$n_sample == 1 && xi_2$n_sample > 1) {
    tmp <- xi_1
    xi_1 <- xi_2
    xi_2 <- tmp
  }

  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample
  if (!(n1 == n2 || n2 == 1)) {
    stop("Sample sizes are incompatible for broadcasting.", call. = FALSE)
  }

  f1 <- xi_1$functional_list
  f2 <- xi_2$functional_list
  Z1 <- xi_1$Z
  Z2 <- xi_2$Z

  if (n2 == 1 && n1 > 1) {
    f2 <- rep_fd(f2, n1)
    Z2 <- matrix(rep(c(Z2), n1), nrow = n1, byrow = TRUE)
  }

  new_functional_list <- Map(
    function(fd1, fd2) plus.fd(fd1, times.fd(alpha, fd2)),
    f1,
    f2
  )
  new_Z <- Z1 + alpha * Z2

  predictor_hybrid(
    Z = new_Z,
    functional_list = new_functional_list,
    eval_point = xi_1$eval_point
  )
}
