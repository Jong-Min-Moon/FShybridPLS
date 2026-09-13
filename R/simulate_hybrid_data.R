#' Simulate Hybrid Functional-Scalar Regression Data
#'
#' Generates synthetic hybrid predictors (functional curves + scalar covariates)
#' and a scalar response useful for examples and testing.
#'
#' @param n Number of samples.
#' @param n_functional Number of functional predictors.
#' @param n_scalar Number of scalar predictors.
#' @param n_basis Number of B-spline basis functions per functional predictor.
#' @param n_eval Number of evaluation grid points on \eqn{[0,1]}.
#' @param signal_to_noise Approximate signal-to-noise ratio for the response.
#' @param seed Optional RNG seed.
#'
#' @return A list with:
#' \describe{
#'   \item{W}{A `predictor_hybrid` object.}
#'   \item{y}{Numeric response vector of length `n`.}
#'   \item{beta_true}{Single-sample `predictor_hybrid` used to generate the linear signal.}
#' }
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 25, n_functional = 1, n_scalar = 2, n_basis = 5)
#' class(sim$W)
#' length(sim$y)
#'
#' @export
#' @importFrom fda create.bspline.basis fd
#' @importFrom stats rnorm sd
simulate_hybrid_data <- function(n = 50,
                                 n_functional = 1,
                                 n_scalar = 3,
                                 n_basis = 7,
                                 n_eval = 51,
                                 signal_to_noise = 3,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(
    is.numeric(n), n >= 5,
    n_functional >= 1, n_scalar >= 1,
    n_basis >= 4, n_eval >= 10
  )

  eval_point <- seq(0, 1, length.out = n_eval)
  basis <- create.bspline.basis(c(0, 1), nbasis = n_basis)

  functional_list <- vector("list", n_functional)
  for (k in seq_len(n_functional)) {
    coefs <- matrix(rnorm(n_basis * n), n_basis, n)
    functional_list[[k]] <- fd(coefs, basis)
  }
  Z <- matrix(rnorm(n * n_scalar), n, n_scalar)
  colnames(Z) <- paste0("Z", seq_len(n_scalar))

  W <- predictor_hybrid(
    Z = Z,
    functional_list = functional_list,
    eval_point = eval_point
  )

  # True coefficient direction in the same hybrid space
  beta_coefs <- c(
    as.numeric(matrix(rnorm(n_basis * n_functional), n_basis, n_functional)),
    rnorm(n_scalar)
  )
  beta_true <- predictor_hybrid_from_coef(format = W, coef = beta_coefs)

  signal <- inprod_predictor_hybrid(W, beta_true)
  noise_sd <- sd(signal) / signal_to_noise
  if (!is.finite(noise_sd) || noise_sd == 0) noise_sd <- 1
  y <- as.numeric(signal + rnorm(n, sd = noise_sd))

  list(W = W, y = y, beta_true = beta_true)
}