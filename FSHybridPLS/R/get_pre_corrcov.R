# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Compute V matrix (pre-corrected covariance estimator)
#'
#' Constructs the V matrix involving inner products of functional predictors
#' (weighted by their Gram matrices) and scalar predictors against response y.
#'
#' @param W A `predictor_hybrid` object.
#' @param y A numeric response vector of length equal to number of samples.
#'
#' @return A matrix representing the corrected cross-covariance structure.
#' @export
get_pre_corrcov <- function(W, y){
  n <- W$n_sample
  K <- W$n_functional

  # building blocks
  ## C_J : C multiplied by J. Scalable to arbitrary K
  predictor_first <- W$functional_list[[1]]
  C_J <- t(predictor_first$coefs) %*% W$gram_list[[1]]
  for (i in 2:K){
    predictor_now <- W$functional_list[[i]]
    C_J <- cbind(
      C_J,
      t(predictor_now$coefs) %*% W$gram_list[[i]]
    )
    }

  
  ## covariance
  J_Ct_y <- t(C_J) %*% y
  Zt_y <- t(W$Z) %*% y

  # V matrix
  upper_left <- J_Ct_y %*% t(J_Ct_y)
  upper_right <- J_Ct_y %*% t(Zt_y)
  lower_left <- t(upper_right)
  lower_right <- Zt_y %*% t(Zt_y)

  V_star <- cbind(
    rbind(upper_left, lower_left),
    rbind(upper_right, lower_right)
    )/(n^2)

  return(V_star)
  }
