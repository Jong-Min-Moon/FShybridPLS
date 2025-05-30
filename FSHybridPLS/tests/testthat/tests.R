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

testthat::test_that("inprod.predictor_hybrid works with broadcasting and self-inner product", {
  suppressPackageStartupMessages(library(fda))
  
  basis <- create.bspline.basis(c(0, 1), 5)
  coef_val <- function(val) matrix(val, 5, 3)

  # Define functional objects
  fd_1 <- fd(coef = coef_val(1), basisobj = basis)
  fd_2 <- fd(coef = coef_val(2), basisobj = basis)

  # Define scalar matrices
  Z1 <- matrix(c(1, 2), nrow = 1)         # 1 × 2
  Z2 <- matrix(c(3, 4), nrow = 1)         # 1 × 2
  Z_many <- matrix(rep(c(1, 2), 3), nrow = 3, byrow = TRUE)  # 3 × 2

  # (1, 1)
  x1 <- predictor_hybrid(Z1, list(fd_1, fd_1), list(matrix(0), matrix(0)), c(5, 5), 1, 2, 2)
  x2 <- predictor_hybrid(Z2, list(fd_2, fd_2), list(matrix(0), matrix(0)), c(5, 5), 1, 2, 2)
  val_scalar <- inprod.predictor_hybrid(x1, x2)

  testthat::expect_type(val_scalar, "double")
  testthat::expect_length(val_scalar, 1)

  # Expected functional inner product: 2 fds × sum(1 * 2) = 2 × (5×3) = 30
  expected_functional <- 2 * sum(1 * 2) * 15  # 5 × 3 = 15 elements
  expected_scalar <- sum(Z1 * Z2)            # 1*3 + 2*4 = 11
  testthat::expect_equal(val_scalar, expected_functional + expected_scalar)

  # (n, 1) broadcasting
  x3 <- predictor_hybrid(Z_many, list(fd_1, fd_1), list(matrix(0), matrix(0)), c(5, 5), 3, 2, 2)
  mat_3_1 <- inprod.predictor_hybrid(x3, x2)
  testthat::expect_equal(dim(mat_3_1), c(3, 1))
  testthat::expect_equal(as.numeric(mat_3_1), rep(val_scalar, 3))

  # (1, m) broadcasting
  x4 <- predictor_hybrid(Z_many, list(fd_2, fd_2), list(matrix(0), matrix(0)), c(5, 5), 3, 2, 2)
  mat_1_3 <- inprod.predictor_hybrid(x1, x4)
  testthat::expect_equal(dim(mat_1_3), c(1, 3))
  testthat::expect_equal(as.numeric(mat_1_3), rep(val_scalar, 3))

  # (n, m)
  mat_3_3 <- inprod.predictor_hybrid(x3, x4)
  testthat::expect_equal(dim(mat_3_3), c(3, 3))
  testthat::expect_true(all(as.numeric(mat_3_3) == val_scalar))

  # shorthand inprod(x) ≡ inprod(x, x)
  self_inner <- inprod.predictor_hybrid(x1)
  testthat::expect_equal(
    self_inner,
    inprod.predictor_hybrid(x1, x1)
  )
})


