#' FSHybridPLS: Hybrid Penalized Partial Least Squares
#'
#' @description
#' Implements Hybrid Penalized Partial Least Squares regression for predictors
#' that combine functional (`fd`) curves and scalar covariates in a joint
#' Hilbert space, with roughness penalties on functional coefficient directions.
#'
#' The main workflow is:
#' 1. Build a hybrid predictor with [predictor_hybrid()] (or
#'    [simulate_hybrid_data()] for synthetic examples).
#' 2. Optionally preprocess with [split_and_normalize_all()].
#' 3. Fit with [fit_hybridPLS()] (optionally tune components via
#'    [cv_fit_hybridPLS()]).
#' 4. Predict with [predict.hybridPLS()].
#'
#' The method is described in Mun and Jang (2026)
#' <doi:10.48550/arXiv.2601.16364>.
#'
#' @keywords internal
"_PACKAGE"
