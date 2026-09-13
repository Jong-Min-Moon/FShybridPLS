#' Compute Penalized PLS Weight Vector (Linear)
#'
#' Solves for the PLS weight direction xi that maximizes the covariance with y,
#' subject to roughness penalties on the functional components.
#'
#' @param W A `predictor_hybrid` object.
#' @param y A numeric vector of response values.
#' @param lambda A numeric vector of smoothing parameters.
#'
#' @return A single-sample `predictor_hybrid` object representing the weight direction xi.
#' @keywords internal
get_xi_hat_linear_pen <- function(W, y, lambda) {
  K <- W$n_functional
  if (length(lambda) != K) {
    stop("'lambda' must have length W$n_functional.", call. = FALSE)
  }
  if (any(diff(W$n_basis_list) != 0)) {
    stop("All functional predictors must share the same number of basis functions.", call. = FALSE)
  }

  u <- gamma <- vector("list", K)
  v <- as.numeric(t(W$Z) %*% y)
  q <- sum(v^2)

  for (j in seq_len(K)) {
    Theta_t <- W$functional_list[[j]]$coefs
    B <- W$gram_list[[j]]
    u[[j]] <- B %*% Theta_t %*% y
    R <- W$gram_list[[j]] + lambda[j] * W$gram_deriv_list[[j]]
    gamma[[j]] <- tryCatch(
      solve(R, u[[j]]),
      error = function(e) {
        stop(
          sprintf("Penalized Gram system for functional predictor %d is singular.", j),
          call. = FALSE
        )
      }
    )
    q <- q + sum(gamma[[j]] * (R %*% gamma[[j]]))
  }

  if (!is.finite(q) || q <= .Machine$double.eps) {
    stop("Penalized weight norm is zero; cannot normalize xi.", call. = FALSE)
  }

  d_vec <- c(do.call(c, gamma), v)
  d_vec <- d_vec / sqrt(q)
  predictor_hybrid_from_coef(format = W, coef = d_vec)
}
