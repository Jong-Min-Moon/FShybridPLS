#' Load and Preprocess Kidney Data for FSHybridPLS
#'
#' Processes a raw kidney study data frame into a `predictor_hybrid` object and
#' a scaled response. The input is expected to contain at least columns
#' `ID`, `Study`, and `renogram_value`, plus clinical covariates used for the
#' response and scalar predictors (accessed by the positional columns described
#' in the paper workflow).
#'
#' This helper is application-specific and requires user-supplied data; it is
#' not needed for the core Hybrid PLS algorithm.
#'
#' @param kidney_df Data frame with kidney study records.
#' @param n_basis Integer. Number of B-spline basis functions (default 20).
#'
#' @return A list with:
#' \describe{
#'   \item{W}{A `predictor_hybrid` object.}
#'   \item{y}{Min-max scaled numeric response in \eqn{[0,1]}.}
#' }
#'
#' @examples
#' try(load_and_preprocess_kidney_data(data.frame(x = 1)), silent = TRUE)
#'
#' @export
#' @importFrom fda create.bspline.basis Data2fd
load_and_preprocess_kidney_data <- function(kidney_df, n_basis = 20) {
  if (!is.data.frame(kidney_df)) {
    stop("'kidney_df' must be a data frame.", call. = FALSE)
  }
  required <- c("ID", "Study", "renogram_value")
  missing_cols <- setdiff(required, names(kidney_df))
  if (length(missing_cols)) {
    stop(
      sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }
  if (ncol(kidney_df) < 29) {
    stop("'kidney_df' must have at least 29 columns for the paper workflow.", call. = FALSE)
  }
  if (!is.numeric(n_basis) || n_basis < 4) {
    stop("'n_basis' must be an integer >= 4.", call. = FALSE)
  }

  ids <- unique(kidney_df$ID)
  n_sample <- length(ids)

  base_data_idx <- which(kidney_df$Study == "Baseline")
  if (!length(base_data_idx)) {
    stop("No rows with Study == 'Baseline' found.", call. = FALSE)
  }
  patient_meta <- kidney_df[base_data_idx, ]
  patient_meta <- patient_meta[!duplicated(patient_meta$ID), ]

  y_mat <- as.matrix(patient_meta[, c(13, 14, 15)])
  y_raw <- rowMeans(y_mat, na.rm = TRUE)
  y_min <- min(y_raw, na.rm = TRUE)
  y_max <- max(y_raw, na.rm = TRUE)
  if (y_max == y_min) {
    warning("Response variable has zero variance. Skipping min-max scaling.")
    y <- y_raw
  } else {
    y <- (y_raw - y_min) / (y_max - y_min)
  }

  Z <- as.matrix(patient_meta[, c(12, 16:29)])
  colnames(Z) <- c("Age", colnames(patient_meta)[16:29])

  df_base <- kidney_df[kidney_df$Study == "Baseline", ]
  df_post <- kidney_df[kidney_df$Study != "Baseline", ]

  reno_base_mat <- matrix(df_base$renogram_value, nrow = n_sample, byrow = TRUE)
  reno_post_mat <- matrix(df_post$renogram_value, nrow = n_sample, byrow = TRUE)

  max_base <- apply(reno_base_mat, 1, max)
  max_base[max_base == 0] <- 1
  reno_base_norm <- reno_base_mat / max_base
  reno_post_norm <- reno_post_mat / max_base

  t_base <- seq(0, 1, length.out = ncol(reno_base_norm))
  t_post <- seq(0, 1, length.out = ncol(reno_post_norm))

  basis_base <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)
  basis_post <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)

  fd_base <- Data2fd(t_base, t(reno_base_norm), basis_base)
  fd_post <- Data2fd(t_post, t(reno_post_norm), basis_post)
  fd_base$fdnames <- list("Time", "Patient", "Baseline Renogram")
  fd_post$fdnames <- list("Time", "Patient", "Post-Furosemide Renogram")

  W <- predictor_hybrid(
    Z = Z,
    functional_list = list(fd_base, fd_post),
    eval_point = seq(0, 1, length.out = 100)
  )

  list(W = W, y = as.numeric(y))
}
