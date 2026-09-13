## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## ----simulate-----------------------------------------------------------------
library(FSHybridPLS)
set.seed(1)

sim <- simulate_hybrid_data(
  n = 50,
  n_functional = 1,
  n_scalar = 3,
  n_basis = 6
)

sim$W
length(sim$y)

## ----preprocess---------------------------------------------------------------
prep <- split_and_normalize_all(sim$W, sim$y, train_ratio = 0.7)
prep$predictor_train$n_sample
prep$predictor_test$n_sample

## ----fit----------------------------------------------------------------------
fit <- fit_hybridPLS(
  prep$predictor_train,
  prep$response_train,
  n_iter = 3,
  lambda = 1e-3,
  validation_data = list(
    W_test = prep$predictor_test,
    y_test = prep$response_test
  )
)

fit
preds <- predict(fit, prep$predictor_test)
rmse <- sqrt(mean((prep$response_test - preds)^2))
rmse

## ----cv-----------------------------------------------------------------------
cv <- cv_fit_hybridPLS(
  prep$predictor_train,
  prep$response_train,
  n_iter = 4,
  lambda = 1e-3,
  n_fold = 3,
  seed = 1
)
cv$rmse_by_component
cv$best_n_iter

## ----session------------------------------------------------------------------
sessionInfo()

