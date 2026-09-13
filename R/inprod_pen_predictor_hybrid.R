#' Penalized Inner Product for Hybrid Objects
#'
#' Computes the unpenalized hybrid inner product plus roughness penalties
#' on functional components: \eqn{\langle \xi_1, \xi_2 \rangle +
#' \sum_k \lambda_k \langle D^m \xi_{1k}, D^m \xi_{2k} \rangle}.
#'
#' @param xi_1 A `predictor_hybrid` object.
#' @param xi_2 Another `predictor_hybrid` object (optional; defaults to `xi_1`).
#' @param lambda Non-negative numeric vector of length `xi_1$n_functional`.
#'
#' @return A numeric scalar (or vector under broadcasting), the penalized
#'   inner product.
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 20, n_basis = 5)
#' fit <- fit_hybridPLS(sim$W, sim$y, n_iter = 2, lambda = 1e-3)
#' inprod_pen_predictor_hybrid(fit$xi[[1]], fit$xi[[1]], lambda = 1e-3)
#'
#' @export
inprod_pen_predictor_hybrid <- function(xi_1, xi_2 = NULL, lambda) {
  assert_predictor_hybrid(xi_1, "xi_1")
  if (is.null(xi_2)) xi_2 <- xi_1
  assert_predictor_hybrid(xi_2, "xi_2")
  if (!is.numeric(lambda) || length(lambda) != xi_1$n_functional || any(lambda < 0)) {
    stop("'lambda' must be a non-negative vector of length xi_1$n_functional.", call. = FALSE)
  }

  base_val <- inprod_predictor_hybrid(xi_1, xi_2)

  # Broadcasting alignment (same rules as inprod_predictor_hybrid)
  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample
  if (n1 == 1 && n2 > 1) {
    tmp <- xi_1
    xi_1 <- xi_2
    xi_2 <- tmp
    n1 <- xi_1$n_sample
    n2 <- xi_2$n_sample
  }
  is_broadcasting <- (n2 == 1 && n1 > 1)

  pen_term <- numeric(n1)
  for (k in seq_len(xi_1$n_functional)) {
    R <- xi_1$gram_deriv_list[[k]]
    C1 <- xi_1$functional_list[[k]]$coefs
    C2 <- xi_2$functional_list[[k]]$coefs
    R_C2 <- R %*% C2
    if (is_broadcasting) {
      pen_term <- pen_term + lambda[k] * as.vector(t(C1) %*% R_C2)
    } else {
      pen_term <- pen_term + lambda[k] * colSums(C1 * R_C2)
    }
  }

  base_val + pen_term
}
