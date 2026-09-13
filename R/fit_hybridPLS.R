#' Fit Hybrid Penalized PLS
#'
#' Runs the iterative Hybrid Penalized Partial Least Squares algorithm.
#' If `validation_data` is provided, RMSE on the validation set is stored
#' for each extracted component.
#'
#' @param W A `predictor_hybrid` object (training predictors).
#' @param y Numeric vector of training responses (length `W$n_sample`).
#' @param n_iter Integer number of PLS components to extract (>= 1).
#' @param lambda Numeric vector of smoothing parameters, one per functional
#'   predictor (`length(lambda) == W$n_functional`).
#' @param validation_data Optional list with elements `W_test`
#'   (`predictor_hybrid`) and `y_test` (numeric vector).
#'
#' @return An object of class `hybridPLS`, a list with:
#' \describe{
#'   \item{rho}{List of score vectors (length `n_iter`).}
#'   \item{xi}{List of weight directions as single-sample `predictor_hybrid` objects.}
#'   \item{beta}{List of cumulative regression directions; `beta[[L]]` uses the
#'     first `L` components.}
#'   \item{W}{List of deflated predictors; `W[[1]]` is the input, `W[[l+1]]`
#'     is after extracting component `l`.}
#'   \item{delta}{List of predictor loadings.}
#'   \item{nu}{List of response loadings (scalars).}
#'   \item{n_iter, lambda}{Fitting settings.}
#'   \item{validation_rmse}{Numeric vector of validation RMSE by component, or `NULL`.}
#' }
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 40, n_functional = 1, n_scalar = 2, n_basis = 5)
#' fit <- fit_hybridPLS(sim$W, sim$y, n_iter = 2, lambda = 1e-3)
#' fit
#' pred <- predict(fit, sim$W, n_components = 2)
#' cor(pred, sim$y)
#'
#' @export
fit_hybridPLS <- function(W, y, n_iter, lambda, validation_data = NULL) {
  assert_predictor_hybrid(W, "W")
  if (!is.numeric(y) || length(y) != W$n_sample) {
    stop("'y' must be a numeric vector of length W$n_sample.", call. = FALSE)
  }
  if (!is.numeric(n_iter) || length(n_iter) != 1 || n_iter < 1 || n_iter != as.integer(n_iter)) {
    stop("'n_iter' must be a positive integer.", call. = FALSE)
  }
  n_iter <- as.integer(n_iter)
  if (!is.numeric(lambda) || length(lambda) != W$n_functional || any(lambda < 0)) {
    stop("'lambda' must be a non-negative numeric vector of length W$n_functional.", call. = FALSE)
  }
  if (!is.null(validation_data)) {
    if (!is.list(validation_data) ||
        is.null(validation_data$W_test) ||
        is.null(validation_data$y_test)) {
      stop("'validation_data' must be a list with W_test and y_test.", call. = FALSE)
    }
    assert_predictor_hybrid(validation_data$W_test, "validation_data$W_test")
    if (!is.numeric(validation_data$y_test) ||
        length(validation_data$y_test) != validation_data$W_test$n_sample) {
      stop("'validation_data$y_test' length must match W_test$n_sample.", call. = FALSE)
    }
  }

  W_now <- rho <- xi <- delta <- nu <- iota <- beta <- vector("list", n_iter + 1)
  validation_rmse <- if (!is.null(validation_data)) numeric(n_iter) else NULL

  W_now[[1]] <- W
  y_now <- y

  for (l in seq_len(n_iter)) {
    xi[[l]] <- get_xi_hat_linear_pen(W_now[[l]], y_now, lambda)
    rho[[l]] <- inprod_predictor_hybrid(W_now[[l]], xi[[l]])
    delta[[l]] <- get_delta(W_now[[l]], rho[[l]])
    nu[[l]] <- get_nu(y_now, rho[[l]])

    W_now[[l + 1]] <- residualize_predictor(W_now[[l]], rho[[l]], delta[[l]])
    y_now <- residualize_y(y_now, rho[[l]], nu[[l]])

    iota[[l]] <- xi[[l]]
    if (l == 1) {
      beta[[l]] <- scalar_mul_predictor_hybrid(iota[[l]], nu[[l]])
    } else {
      for (u in seq_len(l - 1)) {
        adjustment <- inprod_predictor_hybrid(delta[[u]], xi[[l]])
        iota[[l]] <- subtr_predictor_hybrid(iota[[l]], iota[[u]], adjustment)
      }
      beta[[l]] <- add_predictor_hybrid(beta[[l - 1]], iota[[l]], nu[[l]])
    }

    if (!is.null(validation_data)) {
      y_pred_test <- inprod_predictor_hybrid(validation_data$W_test, beta[[l]])
      validation_rmse[l] <- sqrt(mean((validation_data$y_test - y_pred_test)^2))
    }
  }

  structure(
    list(
      rho = rho[seq_len(n_iter)],
      xi = xi[seq_len(n_iter)],
      beta = beta[seq_len(n_iter)],
      W = W_now[seq_len(n_iter + 1)],
      delta = delta[seq_len(n_iter)],
      nu = nu[seq_len(n_iter)],
      n_iter = n_iter,
      lambda = lambda,
      validation_rmse = validation_rmse
    ),
    class = "hybridPLS"
  )
}

#' @param x A `hybridPLS` object.
#' @param ... Unused.
#' @rdname fit_hybridPLS
#' @export
#' @method print hybridPLS
print.hybridPLS <- function(x, ...) {
  cat("Hybrid Penalized PLS model\n")
  cat(sprintf("  Components (n_iter): %d\n", x$n_iter))
  cat(sprintf("  Functional predictors: %d\n", length(x$lambda)))
  cat(sprintf("  Lambda: %s\n", paste(signif(x$lambda, 4), collapse = ", ")))
  if (!is.null(x$validation_rmse)) {
    cat("  Validation RMSE by component:\n")
    cat(sprintf("    %s\n", paste(signif(x$validation_rmse, 4), collapse = ", ")))
  }
  invisible(x)
}

#' Predict from a Hybrid PLS Model
#'
#' Computes fitted or predicted responses using the cumulative coefficient
#' direction `beta[[n_components]]`.
#'
#' @param object A `hybridPLS` object from [fit_hybridPLS()].
#' @param newdata A `predictor_hybrid` object. If omitted, predictions use the
#'   training predictors stored as `object$W[[1]]`.
#' @param n_components Number of components to use (default: all).
#' @param ... Unused.
#'
#' @return Numeric vector of predictions (length `newdata$n_sample`).
#'
#' @examples
#' set.seed(1)
#' sim <- simulate_hybrid_data(n = 30, n_basis = 5)
#' fit <- fit_hybridPLS(sim$W, sim$y, n_iter = 2, lambda = 1e-3)
#' head(predict(fit, sim$W))
#'
#' @export
#' @method predict hybridPLS
predict.hybridPLS <- function(object, newdata = NULL, n_components = object$n_iter, ...) {
  if (!inherits(object, "hybridPLS")) {
    stop("'object' must inherit from class 'hybridPLS'.", call. = FALSE)
  }
  if (is.null(newdata)) {
    newdata <- object$W[[1]]
  }
  assert_predictor_hybrid(newdata, "newdata")
  if (!is.numeric(n_components) || length(n_components) != 1 ||
      n_components < 1 || n_components > object$n_iter) {
    stop("'n_components' must be between 1 and object$n_iter.", call. = FALSE)
  }
  n_components <- as.integer(n_components)
  inprod_predictor_hybrid(newdata, object$beta[[n_components]])
}