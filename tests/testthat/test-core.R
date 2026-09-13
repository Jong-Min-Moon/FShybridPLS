test_that("simulate and predictor_hybrid construct valid objects", {
  sim <- simulate_hybrid_data(n = 20, n_functional = 1, n_scalar = 2, n_basis = 5, seed = 1)
  expect_s3_class(sim$W, "predictor_hybrid")
  expect_equal(sim$W$n_sample, 20)
  expect_equal(sim$W$n_scalar, 2)
  expect_equal(length(sim$y), 20)
  expect_equal(length(sim$W$gram_list), 1)
})

test_that("fit_hybridPLS returns hybridPLS and predict works", {
  sim <- simulate_hybrid_data(n = 30, n_basis = 5, seed = 2)
  fit <- fit_hybridPLS(sim$W, sim$y, n_iter = 2, lambda = 1e-3)
  expect_s3_class(fit, "hybridPLS")
  expect_equal(fit$n_iter, 2)
  expect_length(fit$beta, 2)
  pred <- predict(fit, sim$W, n_components = 2)
  expect_equal(length(pred), 30)
  expect_true(is.numeric(pred))
  expect_true(cor(pred, sim$y) > 0.2)
})

test_that("split_and_normalize_all produces nonempty splits", {
  sim <- simulate_hybrid_data(n = 40, n_basis = 5, seed = 3)
  prep <- split_and_normalize_all(sim$W, sim$y, train_ratio = 0.7)
  expect_true(prep$predictor_train$n_sample >= 1)
  expect_true(prep$predictor_test$n_sample >= 1)
  expect_equal(
    prep$predictor_train$n_sample + prep$predictor_test$n_sample,
    40
  )
})

test_that("cv_fit_hybridPLS returns best_n_iter", {
  sim <- simulate_hybrid_data(n = 36, n_basis = 5, seed = 4)
  cv <- cv_fit_hybridPLS(sim$W, sim$y, n_iter = 3, lambda = 1e-3, n_fold = 3, seed = 4)
  expect_true(cv$best_n_iter %in% 1:3)
  expect_equal(length(cv$rmse_by_component), 3)
})

test_that("fit_hybridPLS validates inputs", {
  sim <- simulate_hybrid_data(n = 20, n_basis = 5, seed = 5)
  expect_error(fit_hybridPLS(sim$W, sim$y, n_iter = 2, lambda = c(1, 2)), "lambda")
  expect_error(fit_hybridPLS(sim$W, sim$y[1:5], n_iter = 1, lambda = 1e-3), "y")
})
