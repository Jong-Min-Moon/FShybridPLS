# ==============================================================================
# 6_1_geometric_validation.R
# Section 6.1
# Validates Geometric Properties for ONE replication.
# ==============================================================================
# Installation of FSHybridPLS package
# This package is generated using litr, a literate programming tool developed by Jacob Bien and Patrick Vossler (https://jacobbien.github.io/litr-project/). I used litr to build the package and you don't have to install litr just to use the package. To install the package:
#
# devtools::install("your_downloaded_directory/FSHybridPLS")

library(fda)
library(dplyr)
library(tidyr)
library(FSHybridPLS)

# ------------------------------------------------------------------------------
# 0. PARSE ARGUMENTS & SETUP
# ------------------------------------------------------------------------------
# Default User Configuration (Editable)
REPO_DIR <- "."                     # Directory where the script is run (usually root of repo or numerical_studies)
SAVE_DIR_REL <- "./results/geometric" # Relative path for results
REP_ID_DEFAULT <- 1                 # Default replication ID

# Create 'results' directory if it doesn't exist
if (!dir.exists(file.path(REPO_DIR, "results"))) dir.create(file.path(REPO_DIR, "results"))
save_dir <- file.path(REPO_DIR, SAVE_DIR_REL)
if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# Parse Command Line Arguments (Optional, overrides defaults)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 1) {
  rep_id <- as.integer(args[1])
  message(sprintf("Using Rep ID from Command Line: %d", rep_id))
} else {
  rep_id <- REP_ID_DEFAULT
  message(sprintf("Using Default Rep ID: %d (Run with 'Rscript script.R <RepID>' to change)", rep_id))
}

# Configuration
n_comp <- 6

# Set seed for reproducibility specific to this job ID
set.seed(rep_id * 2025)

message(sprintf("=== Starting Geometric Validation | Rep: %d ===", rep_id))

# ------------------------------------------------------------------------------
# 1. DATA GENERATION (Latent Factor Model)
# ------------------------------------------------------------------------------
generate_torture_data <- function(n_sample = 100) {
  # 1. Latent Factors (Scale Mismatch 100:1)
  t1 <- rnorm(n_sample, sd = 10)   # High variance
  t2 <- rnorm(n_sample, sd = 0.1)  # Low variance

  # 2. Functional Predictors (Basis M=15)
  eval_points <- seq(0, 1, length.out = 100)
  basis <- create.bspline.basis(c(0, 1), 15)

  # X1: Pure signal from t1
  X1_mat <- outer(t1, sin(2 * pi * eval_points))
  fd1 <- Data2fd(eval_points, t(X1_mat), basis)

  # X2: Signal from t2 + Noise
  noise_func <- matrix(rnorm(n_sample * 100, sd = 0.01), 100, n_sample)
  X2_mat <- outer(t2, sin(10 * pi * eval_points)) + noise_func
  fd2 <- Data2fd(eval_points, t(X2_mat), basis)

  # 3. Rank-Deficient Scalar Predictors (Z)
  P <- 50
  A <- matrix(runif(P * 2, min = -1, max = 1), nrow = P, ncol = 2)
  T_factors <- cbind(t1, t2)
  Z_clean <- T_factors %*% t(A)
  Z_noise <- matrix(rnorm(n_sample * P, sd = 1), n_sample, P)
  Z <- Z_clean + Z_noise

  # 4. Response Y
  y_clean <- 0.5 * t1 + 10 * t2
  y <- y_clean + rnorm(n_sample, sd = 1)

  W_init <- predictor_hybrid(Z, list(fd1, fd2), eval_points)
  return(list(W = W_init, y = y))
}

# ------------------------------------------------------------------------------
# 2. SCENARIOS
# ------------------------------------------------------------------------------
scenarios <- list(
  "Weak"   = c(0.1, 0.1),
  "Mixed"  = c(0.1, 10.0),
  "Strong" = c(10.0, 10.0)
)

# Lists to collect results for this single rep
list_ortho_weights <- list()
list_ortho_scores  <- list()
list_norms         <- list()
list_y_corrs <- list()
# ------------------------------------------------------------------------------
# 3. MAIN LOOP (Through Scenarios only)
# ------------------------------------------------------------------------------

for(scen_name in names(scenarios)) {
  lam_vec <- scenarios[[scen_name]]


  data <- generate_torture_data(100)
  processed <- split_and_normalize.all(data$W, data$y, train_ratio=0.99)
  W_train <- processed$predictor_train
  y_train <- processed$response_train

  # B. Fit Model
  fit <- fit.hybridPLS(W_train, y_train, n_iter=n_comp, lambda=lam_vec)

  # --- ANALYSIS 1: WEIGHT ORTHONORMALITY ---
  for(r in 1:n_comp) {
    for(c in 1:n_comp) {
      val_pen <- inprod_pen.predictor_hybrid(fit$xi[[r]], fit$xi[[c]], lam_vec)
      val_unpen <- inprod.predictor_hybrid(fit$xi[[r]], fit$xi[[c]])

      list_ortho_weights[[length(list_ortho_weights) + 1]] <- data.frame(
        Scenario = scen_name,
        Rep = rep_id,
        Row = r,
        Col = c,
        Type = "Weight",
        Penalized_InnerProd = val_pen,
        Unpenalized_InnerProd = val_unpen,
        Is_Diagonal = (r == c)
      )
    }
  }

  # --- ANALYSIS 2: SCORE ORTHOGONALITY ---
  for(r in 1:n_comp) {
    for(c in 1:n_comp) {
      val_score_dot <- sum(fit$rho[[r]] * fit$rho[[c]])
      val_score_cor <- cor(fit$rho[[r]], fit$rho[[c]])
      if(is.na(val_score_cor)) val_score_cor <- 0

      list_ortho_scores[[length(list_ortho_scores) + 1]] <- data.frame(
        Scenario = scen_name,
        Rep = rep_id,
        Row = r,
        Col = c,
        Type = "Score",
        Dot_Product = val_score_dot,
        Correlation = val_score_cor,
        Is_Diagonal = (r == c)
      )
    }
  }

  # --- ANALYSIS 3: CROSS-PROJECTION NORMS  --
  for(l in 1:(n_comp-1)) {
    for(k in (l+1):n_comp) {
      proj_vec <- inprod.predictor_hybrid(fit$W[[k]], fit$xi[[l]])
      l2_norm_val <- sqrt(sum(proj_vec^2))

      list_norms[[length(list_norms) + 1]] <- data.frame(
        Scenario = scen_name,
        Rep = rep_id,
        Predictor_Step_L = l,
        Weight_Step_K = k,
        Sample_L2_Norm = l2_norm_val
      )
    }
  }
  # ... inside the loop, after fitting the model ...

    # --- ANALYSIS 4: Y-SCORE CORRELATION ---
    Rho_mat <- do.call(cbind, fit$rho)

    # Calculate correlation between Y (training) and each Score vector
    # We expect high correlation for early components
    y_corrs <- cor(y_train, Rho_mat)

    for(k in 1:n_comp) {
      # Add to a new list (make sure to initialize list_y_corrs <- list() at start)
      list_y_corrs[[length(list_y_corrs) + 1]] <- data.frame(
        Scenario = scen_name,
        Rep = rep_id,
        Index1 = k,          # Component Index
        Index2 = NA,         # Not applicable
        Measure = "Correlation",
        Category = "Y_Correlations",
        Value = y_corrs[1, k] # Correlation value
      )
    }
}


# ------------------------------------------------------------------------------
# 4. SAVE RESULTS
# ------------------------------------------------------------------------------
df_weights <- bind_rows(list_ortho_weights)
df_scores  <- bind_rows(list_ortho_scores)
df_norms   <- bind_rows(list_norms)
df_y_corrs <- bind_rows(list_y_corrs)  # <--- Added this

output_filename <- file.path(save_dir, sprintf("geom_results_rep_%d.rds", rep_id))

saveRDS(list(
  weights = df_weights,
  scores = df_scores,
  norms = df_norms,
  y_corrs = df_y_corrs                 # <--- Added this
), file = output_filename)

message(sprintf("=== Replication %d Complete. Saved to %s ===", rep_id, output_filename))
