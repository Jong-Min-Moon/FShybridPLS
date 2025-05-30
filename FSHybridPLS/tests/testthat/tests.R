# Generated from create-rhello.Rmd: do not edit by hand  
testthat::test_that("predictor_hybrid constructor works as expected (S3)", {
  Z <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  fd1 <- list(fdobj = "fd1 placeholder")
  fd2 <- list(fdobj = "fd2 placeholder")
  J1 <- matrix(runif(100), nrow = 10)
  J2 <- matrix(runif(200), nrow = 20)
  
  obj <- predictor_hybrid(
    Z = Z,
    functional_list = list(fd1, fd2),
    jacobian_list = list(J1, J2),
    n_basis_list = c(10, 20),
    n_sample = 100,
    n_functional = 2,
    n_scalar = 5
  )
  
  testthat::expect_s3_class(obj, "predictor_hybrid")
  testthat::expect_equal(nrow(obj$Z), 100)
  testthat::expect_equal(ncol(obj$Z), 5)
  testthat::expect_equal(length(obj$functional_list), 2)
  testthat::expect_equal(length(obj$jacobian_list), 2)
  testthat::expect_equal(obj$n_basis_list, c(10, 20))
  testthat::expect_equal(obj$n_sample, 100)
  testthat::expect_equal(obj$n_functional, 2)
  testthat::expect_equal(obj$n_scalar, 5)
})



testthat::test_that("add.predictor_hybrid works for matching predictors", {
  suppressPackageStartupMessages(library(fda))

  Z1 <- matrix(1, nrow = 3, ncol = 2)
  Z2 <- matrix(2, nrow = 3, ncol = 2)
  
  basis <- fda::create.bspline.basis(c(0, 1), 5)
  fd1 <- fda::fd(coef = matrix(1, 5, 3), basisobj = basis)
  fd2 <- fda::fd(coef = matrix(2, 5, 3), basisobj = basis)
  
  obj1 <- predictor_hybrid(Z1, list(fd1), list(matrix(0)), c(5), 3, 1, 2)
  obj2 <- predictor_hybrid(Z2, list(fd2), list(matrix(0)), c(5), 3, 1, 2)
  
  obj3 <- add.predictor_hybrid(obj1, obj2, alpha = 1)
  
  testthat::expect_s3_class(obj3, "predictor_hybrid")
  testthat::expect_equal(obj3$Z, matrix(3, nrow = 3, ncol = 2))
  
  expected_coef <- matrix(3, 5, 3)
  testthat::expect_equal(
    unname(coef(obj3$functional_list[[1]])),
    unname(expected_coef)
  )
})



testthat::test_that("scalar_mul.predictor_hybrid works with two functional predictors", {
  suppressPackageStartupMessages(library(fda))
  
  # Scalar matrix
  Z <- matrix(2, nrow = 3, ncol = 2)
  
  # Create a shared basis for both functional predictors
  basis <- create.bspline.basis(c(0, 1), 5)
  
  # Two identical fd objects
  fd1 <- fd(coef = matrix(2, 5, 3), basisobj = basis)
  fd2 <- fd(coef = matrix(2, 5, 3), basisobj = basis)
  
  # Construct predictor_hybrid object
  obj <- predictor_hybrid(
    Z = Z,
    functional_list = list(fd1, fd2),
    jacobian_list = list(matrix(0), matrix(0)),
    n_basis_list = c(5, 5),
    n_sample = 3,
    n_functional = 2,
    n_scalar = 2
  )
  
  # Perform scalar multiplication
  scaled <- scalar_mul.predictor_hybrid(obj, 4)
  
  testthat::expect_s3_class(scaled, "predictor_hybrid")
  
  # Scalar part should be scaled
  testthat::expect_equal(scaled$Z, Z * 4)
  
  # Each functional part should be scaled
  expected_coef <- matrix(8, 5, 3)  # 2 * 4
  for (fd_scaled in scaled$functional_list) {
    testthat::expect_equal(
      unname(coef(fd_scaled)),
      unname(expected_coef)
    )
  }
})

