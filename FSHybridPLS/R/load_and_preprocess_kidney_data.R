# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Load and Preprocess Kidney Data for FSHybridPLS
#'
#' This function processes the raw 'kidney' dataframe into a `predictor_hybrid` object.
#' It performs:
#' 1. Extraction of scalar and functional covariates.
#' 2. Specific medical normalization (scaling by baseline max).
#' 3. Smoothing of functional curves into B-spline basis objects.
#' 4. Transformation of Response: Min-Max scaling followed by Logit transformation.
#' 5. Construction of the hybrid predictor object.
#'
#' @param kidney_df The raw dataframe containing kidney data (columns: ID, Study, renogram_value, etc.)
#' @param n_basis Integer. Number of B-spline basis functions to use for smoothing (default 20).
#'
#' @return A list containing:
#'   /item{W}{A `predictor_hybrid` object containing the predictors.}
#'   /item{y}{A numeric vector representing the transformed response (logit of min-max scaled mean diagnosis metrics).}
#' @export
load_and_preprocess_kidney_data <- function(kidney_df, n_basis = 20) {
  
  # --- 1. Data Extraction ---
  # Split data by ID to ensure order consistency
  ids <- unique(kidney_df$ID)
  n_sample <- length(ids)
  
  # Helper to extract scalar/response from the first row of each patient's baseline
  base_data_idx <- which(kidney_df$Study == "Baseline")
  patient_meta <- kidney_df[base_data_idx, ]
  patient_meta <- patient_meta[!duplicated(patient_meta$ID), ]
  
  # --- Response Extraction & Transformation ---
  # 1. Raw Mean (Columns 13, 14, 15)
  y_mat <- as.matrix(patient_meta[, c(13, 14, 15)])
  y_raw <- rowMeans(y_mat, na.rm = TRUE)
  
  # 2. Min-Max Scaling to [0, 1]
  y_min <- min(y_raw, na.rm = TRUE)
  y_max <- max(y_raw, na.rm = TRUE)
  # Prevent division by zero if all y are identical
  if (y_max == y_min) {
    warning("Response variable has zero variance. Skipping scaling/logit.")
    y <- y_raw
  } else {
    y_scaled <- (y_raw - y_min) / (y_max - y_min)
    
    # 3. Safety Squeeze (avoid Inf/-Inf at boundaries)
    # Maps [0, 1] -> [epsilon, 1-epsilon]
    epsilon <- 1e-5
    y_safe <- y_scaled * (1 - 2 * epsilon) + epsilon
    
    # 4. Logit Transformation: log(p / (1-p))
    #y <- qlogis(y_safe)
    y <- y_scaled
    # <- y_raw
  }
  
  # --- Scalar Predictors (Z) ---
  # Columns: 12 (Age) + 16:29
  Z <- as.matrix(patient_meta[, c(12, 16:29)])
  colnames(Z) <- c("Age", colnames(patient_meta)[16:29])
  
  # --- Functional Data Extraction ---
  # Filter and sort
  df_base <- kidney_df[kidney_df$Study == "Baseline", ]
  df_post <- kidney_df[kidney_df$Study != "Baseline", ]
  
  # Reshape to (ID x Time) matrices
  reno_base_mat <- matrix(df_base$renogram_value, nrow = n_sample, byrow = TRUE)
  reno_post_mat <- matrix(df_post$renogram_value, nrow = n_sample, byrow = TRUE)
  
  # --- 2. Medical Preprocessing (Normalization) ---
  # Calculate max of baseline for each patient
  max_base <- apply(reno_base_mat, 1, max)
  max_base[max_base == 0] <- 1 
  
  # Scale matrices
  reno_base_norm <- reno_base_mat / max_base
  reno_post_norm <- reno_post_mat / max_base
  
  # --- 3. Smoothing (Convert to fd objects) ---
  # Define time grids normalized to [0,1]
  t_base <- seq(0, 1, length.out = ncol(reno_base_norm))
  t_post <- seq(0, 1, length.out = ncol(reno_post_norm))
  
  # Create Basis (B-spline)
  basis_base <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)
  basis_post <- create.bspline.basis(rangeval = c(0, 1), nbasis = n_basis)
  
  # Smooth (Transpose: Time x Samples)
  fd_base <- Data2fd(t_base, t(reno_base_norm), basis_base)
  fd_post <- Data2fd(t_post, t(reno_post_norm), basis_post)
  
  # Assign Names
  fd_base$fdnames <- list("Time", "Patient", "Baseline Renogram")
  fd_post$fdnames <- list("Time", "Patient", "Post-Furosemide Renogram")
  
  # --- 4. Construct Hybrid Predictor ---
  common_eval <- seq(0, 1, length.out = 100)
  
  W <- predictor_hybrid(
    Z = Z,
    functional_list = list(fd_base, fd_post),
    eval_point = common_eval
  )
  
  return(list(W = W, y = y))
}
