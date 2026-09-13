# Basis equality check without touching fda::: internals
is_same_basis <- function(b1, b2) {
  identical(b1$type, b2$type) &&
    isTRUE(all.equal(b1$rangeval, b2$rangeval, tolerance = sqrt(.Machine$double.eps))) &&
    identical(b1$nbasis, b2$nbasis) &&
    isTRUE(all.equal(b1$params, b2$params, tolerance = sqrt(.Machine$double.eps))) &&
    identical(b1$dropind, b2$dropind)
}

assert_predictor_hybrid <- function(x, arg = "W") {
  if (!inherits(x, "predictor_hybrid")) {
    stop(sprintf("'%s' must be a predictor_hybrid object.", arg), call. = FALSE)
  }
  invisible(x)
}