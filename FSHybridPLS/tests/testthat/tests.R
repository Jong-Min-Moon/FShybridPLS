# Generated from create-rhello.Rmd: do not edit by hand  
testthat::test_that("Gram matrix edge cases", {
  suppressPackageStartupMessages(library(fda))

  # Case 1: Fourier basis (orthogonal over [0, 1])
  fb <- create.fourier.basis(rangeval = c(0, 1), nbasis = 5)
  J_fourier <- compute_gram_matrix(fb)
  testthat::expect_equal(J_fourier, diag(diag(J_fourier)), tolerance = 1e-10)  # Nearly diagonal

  # Case 2: Constant basis (1 basis function)
  cb <- create.constant.basis(rangeval = c(0, 1))
  J_const <- compute_gram_matrix(cb)
  testthat::expect_equal(dim(J_const), c(1, 1))
  testthat::expect_equal(J_const[1, 1], 1)

  # Case 3: Invalid input
  testthat::expect_error(compute_gram_matrix("not a basisfd"))

  # Case 4: Highly localized B-splines
  # Simulate nearly disjoint support (narrow B-splines)
  narrow_bs <- create.bspline.basis(rangeval = c(0, 1), nbasis = 10, norder = 2)
  J_narrow <- compute_gram_matrix(narrow_bs)
  testthat::expect_true(all(abs(J_narrow[upper.tri(J_narrow)]) < 0.5))  # off-diagonals small
})



testthat::test_that("rep_fd correctly replicates a list of fd objects", {
  suppressPackageStartupMessages(library(fda))

  # Create a simple basis
  basis <- create.bspline.basis(c(0, 1), 5)

  # Create single-sample fd object
  coef_vec <- 1:5
  fd1 <- fd(coef = matrix(coef_vec, ncol = 1), basisobj = basis)

  # Replicate the list
  rep_list <- rep_fd(list(fd1, fd1), 3)

  # Check: list length unchanged
  testthat::expect_equal(length(rep_list), 2)

  # Check each replicated fd object
  for (fd_obj in rep_list) {
    testthat::expect_s3_class(fd_obj, "fd")
    testthat::expect_equal(dim(fd_obj$coefs), c(5, 3))  # 5 basis × 3 samples
    for (i in 1:3) {
      testthat::expect_equal(fd_obj$coefs[, i], coef_vec)
    }
  }

  # Error if not single-sample
  fd_multi <- fd(coef = matrix(1:10, nrow = 5, ncol = 2), basisobj = basis)
  testthat::expect_error(rep_fd(list(fd_multi)), "Each fd object in the list must have one column")
  
  # Error on non-fd input
  testthat::expect_error(rep_fd(list(1, "not_fd")), "Input must be a list of 'fd' objects")
})



testthat::test_that("predictor_hybrid constructor works as expected (S3, updated)", {
  suppressPackageStartupMessages(library(fda))

  # Create scalar matrix
  Z <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)

  # Create two fd objects with 100 replications
  basis1 <- create.bspline.basis(c(0, 1), nbasis = 10)
  basis2 <- create.bspline.basis(c(0, 1), nbasis = 20)

  coef1 <- matrix(rnorm(10 * 100), nrow = 10)  # 10 basis × 100 samples
  coef2 <- matrix(rnorm(20 * 100), nrow = 20)  # 20 basis × 100 samples

  fd1 <- fd(coef = coef1, basisobj = basis1)
  fd2 <- fd(coef = coef2, basisobj = basis2)

  # Construct object
  obj <- predictor_hybrid(Z = Z, functional_list = list(fd1, fd2))

  # Assertions
  testthat::expect_s3_class(obj, "predictor_hybrid")
  testthat::expect_equal(nrow(obj$Z), 100)
  testthat::expect_equal(ncol(obj$Z), 5)
  testthat::expect_equal(length(obj$functional_list), 2)
  testthat::expect_equal(obj$n_basis_list, c(10, 20))
  testthat::expect_equal(obj$n_sample, 100)
  testthat::expect_equal(obj$n_functional, 2)
  testthat::expect_equal(obj$n_scalar, 5)
  testthat::expect_equal(dim(obj$gram_list[[1]]), c(10, 10))
  testthat::expect_equal(dim(obj$gram_list[[2]]), c(20, 20))
})

testthat::test_that("add.predictor_hybrid handles broadcasting and basis compatibility", {
  suppressPackageStartupMessages(library(fda))

  # Create shared B-spline basis
  basis <- create.bspline.basis(c(0, 1), nbasis = 5)

  # Helper to generate a valid coefficient matrix
  coef_mat <- function(val, n_sample, basis_obj) {
    nbasis <- basis_obj$nbasis
    matrix(val, nrow = nbasis, ncol = n_sample)
  }

  # Helper to create a list of fd objects (2 functional predictors)
  fd_list <- function(val, n_sample, basis_obj = basis) {
    list(
      fd(coef = coef_mat(val, n_sample, basis_obj), basisobj = basis_obj),
      fd(coef = coef_mat(val, n_sample, basis_obj), basisobj = basis_obj)
    )
  }

  # Scalar predictors
  Z1 <- matrix(1, nrow = 1, ncol = 2)
  Z3 <- matrix(3, nrow = 3, ncol = 2)

  # Functional predictors (same basis)
  F1 <- fd_list(1, 1)
  F3 <- fd_list(3, 3)

  # Create predictor_hybrid objects
  x1 <- predictor_hybrid(Z1, F1)  # 1 sample
  x3 <- predictor_hybrid(Z3, F3)  # 3 samples

  # Broadcast x1 to match x3
  result <- add.predictor_hybrid(x1, x3, alpha = 1)

  # Check metadata
  testthat::expect_s3_class(result, "predictor_hybrid")
  testthat::expect_equal(result$n_sample, 3)
  testthat::expect_equal(result$n_scalar, 2)
  testthat::expect_equal(result$n_functional, 2)

  # Check scalar component
  expected_Z <- matrix(4, nrow = 3, ncol = 2)  # 1 + 3 = 4
  testthat::expect_equal(result$Z, expected_Z)

  # Check functional component
  for (fd_out in result$functional_list) {
    testthat::expect_equal(unname(coef(fd_out)), matrix(4, nrow = 5, ncol = 3))
  }

  # Test symmetry: broadcasting x3 + x1
  result2 <- add.predictor_hybrid(x3, x1, alpha = 1)
  testthat::expect_equal(result2$Z, expected_Z)
  for (fd_out in result2$functional_list) {
    testthat::expect_equal(unname(coef(fd_out)), matrix(4, nrow = 5, ncol = 3))
  }

  # Incompatible sample sizes (should fail)
  x_bad <- predictor_hybrid(matrix(1, 2, 2), fd_list(1, 2))
  testthat::expect_error(add.predictor_hybrid(x3, x_bad), "incompatible for broadcasting")

  # Incompatible functional basis (should fail)
  bad_basis <- create.bspline.basis(c(0, 1), nbasis = 6)
  bad_fd <- fd_list(2, 3, bad_basis)
  x_basis_mismatch <- predictor_hybrid(Z3, bad_fd)
  testthat::expect_error(add.predictor_hybrid(x3, x_basis_mismatch), "same basis")
})


testthat::test_that("scalar_mul.predictor_hybrid works with two functional predictors", {
  suppressPackageStartupMessages(library(fda))
  
  # Scalar part (3 observations, 2 scalar predictors)
  Z <- matrix(2, nrow = 3, ncol = 2)

  # Functional part: shared B-spline basis
  basis <- create.bspline.basis(c(0, 1), nbasis = 5)
  fd1 <- fd(coef = matrix(2, nrow = 5, ncol = 3), basisobj = basis)
  fd2 <- fd(coef = matrix(2, nrow = 5, ncol = 3), basisobj = basis)

  # Construct hybrid object using new simplified constructor
  obj <- predictor_hybrid(Z = Z, functional_list = list(fd1, fd2))

  # Apply scalar multiplication
  scaled <- scalar_mul.predictor_hybrid(obj, scalar = 4)

  # Check class
  testthat::expect_s3_class(scaled, "predictor_hybrid")

  # Check scalar predictor was scaled
  testthat::expect_equal(scaled$Z, Z * 4)

  # Check each functional component's coefficients
  expected_coef <- matrix(8, nrow = 5, ncol = 3)
  for (fd_obj in scaled$functional_list) {
    testthat::expect_equal(
      unname(coef(fd_obj)),
      unname(expected_coef)
    )
  }
})

testthat::test_that("inprod.predictor_hybrid returns correct vector output for broadcasting and other cases", {
  suppressPackageStartupMessages(library(fda))

  basis <- create.bspline.basis(c(0, 1), 5)
  make_fd <- function(val, n_sample) {
    fd(coef = matrix(val, 5, n_sample), basisobj = basis)
  }

  # Case 1: Broadcasting (multi vs. single)
  Z1 <- matrix(c(1, 2), nrow = 3, ncol = 2, byrow = TRUE)
  Z2 <- matrix(c(3, 4), nrow = 1)
  fd1 <- make_fd(1, 3)
  fd2 <- make_fd(2, 1)
  x1 <- predictor_hybrid(Z1, list(fd1, fd1))
  x2 <- predictor_hybrid(Z2, list(fd2, fd2))
  out1 <- inprod.predictor_hybrid(x1, x2)
  testthat::expect_type(out1, "double")
  testthat::expect_length(out1, 3)

  # Case 2: Broadcasting (single vs. multi)
  out2 <- inprod.predictor_hybrid(x2, x1)
  testthat::expect_type(out2, "double")
  testthat::expect_length(out2, 3)

  # Case 3: Self inner product (1 sample)
  out3 <- inprod.predictor_hybrid(x2)
  testthat::expect_type(out3, "double")
  testthat::expect_length(out3, 1)

  # Case 4: One-sample vs. one-sample
  Z3 <- matrix(c(5, 6), nrow = 1)
  fd3 <- make_fd(3, 1)
  x3 <- predictor_hybrid(Z3, list(fd3, fd3))
  out4 <- inprod.predictor_hybrid(x2, x3)
  testthat::expect_type(out4, "double")
  testthat::expect_length(out4, 1)

  # Case 5: Two-sample vs. two-sample
  Z4 <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  fd4 <- make_fd(1, 2)
  x4 <- predictor_hybrid(Z4, list(fd4, fd4))
  out5 <- inprod.predictor_hybrid(x4, x4)
  testthat::expect_type(out5, "double")
  testthat::expect_length(out5, 2)
})

