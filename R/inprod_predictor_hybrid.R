#' Algebraic Inner Product for Hybrid Objects
#'
#' Uses pre-computed Gram matrices to calculate hybrid inner products
#' without numerical integration.
#'
#' @param xi_1 A `predictor_hybrid` object.
#' @param xi_2 Another `predictor_hybrid` object (optional). If `NULL`,
#'   the inner product of `xi_1` with itself is computed.
#'
#' @return A numeric scalar or vector of inner products. When one argument
#'   has a single sample and the other has `N` samples, returns a length-`N`
#'   vector (broadcasted projection).
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 20, n_basis = 5)
#' fit <- fit_hybridPLS(sim$W, sim$y, n_iter = 1, lambda = 1e-3)
#' scores <- inprod_predictor_hybrid(sim$W, fit$xi[[1]])
#' length(scores)
#'
#' @export
inprod_predictor_hybrid <- function(xi_1, xi_2 = NULL) {
  assert_predictor_hybrid(xi_1, "xi_1")
  if (is.null(xi_2)) xi_2 <- xi_1
  assert_predictor_hybrid(xi_2, "xi_2")

  n1 <- xi_1$n_sample
  n2 <- xi_2$n_sample
  is_broadcasting <- (n2 == 1 && n1 > 1)

  if (n1 != n2 && !is_broadcasting) {
    if (n1 == 1 && n2 > 1) {
      tmp <- xi_1
      xi_1 <- xi_2
      xi_2 <- tmp
      n1 <- xi_1$n_sample
      n2 <- xi_2$n_sample
      is_broadcasting <- TRUE
    } else {
      stop("Incompatible sample sizes.", call. = FALSE)
    }
  }

  if (is_broadcasting) {
    term_scalar <- as.vector(xi_1$Z %*% t(xi_2$Z))
  } else {
    term_scalar <- rowSums(xi_1$Z * xi_2$Z)
  }

  term_functional <- numeric(n1)
  for (k in seq_len(xi_1$n_functional)) {
    J <- xi_1$gram_list[[k]]
    C1 <- xi_1$functional_list[[k]]$coefs
    C2 <- xi_2$functional_list[[k]]$coefs
    J_C2 <- J %*% C2
    if (is_broadcasting) {
      term_functional <- term_functional + as.vector(t(C1) %*% J_C2)
    } else {
      term_functional <- term_functional + colSums(C1 * J_C2)
    }
  }

  term_scalar + term_functional
}
