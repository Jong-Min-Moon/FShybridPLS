# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Fit Hybrid Penalized PLS with Optional Validation
#'
#' Runs the iterative PLS algorithm. If 'validation_data' is provided,
#' it computes and stores the RMSE on the validation set for each component.
#'
#' @param W Initial `predictor_hybrid` object (Training).
#' @param y Numeric vector of response values (Training).
#' @param n_iter Integer, number of components to extract.
#' @param lambda Numeric vector of smoothing parameters.
#' @param validation_data (Optional) A list containing `W_test` (predictor_hybrid) and `y_test`.
#'
#' @return A list containing the model (rho, xi, beta) and optional `validation_rmse`.
#' @export
fit.hybridPLS <- function(W, y, n_iter, lambda, validation_data = NULL) {
  # 1. Initialize storage
  W_now <- rho <- xi <- delta <- nu <- iota <- beta <- list()
  validation_rmse <- numeric(n_iter)
  
  # 2. Initialize current state
  W_now[[1]] <- W
  y_now <- y
  
  # 3. Main PLS Loop
  for (l in 1:n_iter) {
    # --- Step A: Weight Direction (xi) ---
    xi[[l]] <- get_xi_hat_linear_pen(W_now[[l]], y_now, lambda)
    
    # --- Step B: Scores (rho) ---
    rho[[l]] <- inprod.predictor_hybrid(W_now[[l]], xi[[l]])
    
    # --- Step C: Loadings (delta, nu) ---
    delta[[l]] <- get_delta(W_now[[l]], rho[[l]])
    nu[[l]] <- get_nu(y_now, rho[[l]])
    
    # --- Step D: Deflation ---
    W_now[[l + 1]] <- residualize_predictor(W_now[[l]], rho[[l]], delta[[l]])
    y_now <- residualize_y(y_now, rho[[l]], nu[[l]])
    
    # --- Step E: Cumulative Coefficient (Beta) ---
    iota[[l]] <- xi[[l]]
    
    if (l == 1) {
      beta[[l]] <- scalar_mul.predictor_hybrid(iota[[l]], nu[[l]])
    } else {
      for (u in 1:(l - 1)) {
        adjustment <- inprod.predictor_hybrid(delta[[u]], xi[[l]])
        iota[[l]] <- subtr.predictor_hybrid(iota[[l]], iota[[u]], adjustment)
      }
      beta[[l]] <- add.predictor_hybrid(beta[[l - 1]], iota[[l]], nu[[l]])
    }
    
    # --- Step F: Validation (Optional) ---
    if (!is.null(validation_data)) {
      # Predict on test set using the current cumulative beta
      y_pred_test <- inprod.predictor_hybrid(validation_data$W_test, beta[[l]])
      
      # Compute RMSE (standardized scale usually, but here just raw RMSE)
      validation_rmse[l] <- sqrt(mean((validation_data$y_test - y_pred_test)^2))
    }
  }
  
  return(list(
    rho = rho,
    xi = xi,
    beta = beta,
    validation_rmse = if (!is.null(validation_data)) validation_rmse else NULL
  ))
}
