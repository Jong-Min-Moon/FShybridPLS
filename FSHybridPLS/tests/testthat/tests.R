# Generated from create-FSHybridPLS.Rmd: do not edit by hand  
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

testthat::test_that("hybrid_from_coef reconstructs single-sample predictor_hybrid correctly", {
  suppressPackageStartupMessages(library(fda))

  make_fd <- function(coefs, basis) {
    fd(coef = coefs, basisobj = basis)
  }

  basis <- create.bspline.basis(c(0, 1), nbasis = 4)
  n_sample <- 5
  n_scalar <- 2
  n_functional <- 2

  # Create dummy functional list with zeros
  fd1 <- make_fd(matrix(0, 4, n_sample), basis)
  fd2 <- make_fd(matrix(0, 4, n_sample), basis)
  Z <- matrix(0, n_sample, n_scalar)
  obj_template <- predictor_hybrid(Z = Z, functional_list = list(fd1, fd2))

  # Create coefficient vector for reconstruction: 4 + 4 + 2 = 10
  xi_star <- c(1:10)

  # Apply reconstruction
  obj_single <- predictor_hybrid_from_coef(obj_template, xi_star)

  # Checks
  testthat::expect_s3_class(obj_single, "predictor_hybrid")
  testthat::expect_equal(obj_single$n_sample, 1)
  testthat::expect_equal(obj_single$Z, matrix(c(9, 10), nrow = 1))

  # Check functional coefficients (use coef() not $coefs)
  testthat::expect_equal(drop(coef(obj_single$functional_list[[1]])), 1:4)
  testthat::expect_equal(drop(coef(obj_single$functional_list[[2]])), 5:8)
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
# data generation
## three observations
  fd1 <- make_fd(1, 3)
  Z1 <- matrix(c(1, 2), nrow = 3, ncol = 2, byrow = TRUE)
  x1 <- predictor_hybrid(Z1, list(fd1, fd1))
## one observation
  fd2 <- make_fd(2, 1)
  Z2 <- matrix(c(3, 4), nrow = 1)
  x2 <- predictor_hybrid(Z2, list(fd2, fd2))
## three observations
  Z3 <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
  fd3 <- fd(coef = t(matrix(c(1:15), nrow = 3, byrow = TRUE)), basisobj = basis)
  x3 <- predictor_hybrid(Z3, list(fd3, fd3))  
## one observations
  Z4 <- matrix(c(5, 6), nrow = 1)
  fd4 <- make_fd(3, 1)
  x4 <- predictor_hybrid(Z4, list(fd4, fd4))
  
# Case 1: Self inner product (1 sample)
  out1 <- inprod.predictor_hybrid(x2)
  expected_scalar_1 <- sum(Z2^2)
  expected_functional_1 <- sum(sapply(x2$functional_list, function(fd) fda::inprod(fd, fd)))
  testthat::expect_equal(out1, expected_scalar_1 + expected_functional_1)

# Case 2: Three-sample vs. Three-sample
  out2 <- inprod.predictor_hybrid(x1, x3)
  expected_scalar_2 <- rowSums(Z1 * Z3)
  expected_functional_2 <- vapply(1:3, function(i) {
    sum(sapply(1:2, function(j) {
      fda::inprod(x1$functional_list[[j]][i], x3$functional_list[[j]][i])
    }))
  }, numeric(1))
  testthat::expect_equal(out2, expected_scalar_2 + expected_functional_2)
  
  # Case 3: One-sample vs. one-sample

  out3 <- inprod.predictor_hybrid(x2, x4)
  expected_scalar_3 <- sum(Z2 * Z4)
  expected_functional_3 <-2*fda::inprod(fd2, fd4)
  testthat::expect_equal(out3, as.numeric(expected_scalar_3 + expected_functional_3))


  # Case 4: Broadcasting (single vs. multi)
  out4 <- inprod.predictor_hybrid(x2, x3)
  expected_scalar_4 <- c(0,0,0)
  for (i in 1:3){
    expected_scalar_4[i] <-  sum(Z2 * Z3[i,])
    }

  expected_functional_4 <-  c(0,0,0)
  for (i in 1:3){
    expected_functional_4[i] <- 2*fda::inprod(
      fd2,
      fd(
        coef = matrix(
        c(
        (5*i+1-5):(5*i)
        ), 5, 1), #end of matrix
         basisobj = basis) #end of fd
      )# end of inprod
  }

  testthat::expect_equal(as.vector(out4), as.vector(expected_scalar_4 + expected_functional_4))

  # Case 4: Broadcasting (multi  vs. single)
  testthat::expect_equal(as.vector(inprod.predictor_hybrid( x3, x2)), as.vector(expected_scalar_4 + expected_functional_4))  
})

testthat::test_that("get_gram_matrix_block constructs correct block-diagonal matrix", {
  suppressPackageStartupMessages(library(fda))
  suppressPackageStartupMessages(library(Matrix))

  # Construct hybrid predictor with two functional and three scalar predictors
  basis1 <- create.bspline.basis(c(0, 1), nbasis = 4)
  basis2 <- create.bspline.basis(c(0, 1), nbasis = 4)  # was 3; now valid

  fd1 <- fd(coef = matrix(1, 4, 2), basisobj = basis1)
  fd2 <- fd(coef = matrix(2, 4, 2), basisobj = basis2)
  Z <- matrix(rnorm(2 * 3), nrow = 2, ncol = 3)

  obj <- predictor_hybrid(Z = Z, functional_list = list(fd1, fd2))

  # Extract the block-diagonal Gram matrix
  G <- get_gram_matrix_block(obj)

  # Expected block sizes
  nb1 <- 4
  nb2 <- 4
  ns <- 3
  total_dim <- nb1 + nb2 + ns

  # Check matrix properties
  testthat::expect_equal(dim(G), c(total_dim, total_dim))

  # Check sub-block structure
  # Top-left: Gram matrix 1
  testthat::expect_equal(sum(G[1:nb1, 1:nb1]-obj$gram_list[[1]]), 0)

  # Next block: Gram matrix 2
  testthat::expect_equal(sum( G[(nb1 + 1):(nb1 + nb2), (nb1 + 1):(nb1 + nb2)]-obj$gram_list[[2]]),0 )

  # Final block: identity matrix for scalar part
  testthat::expect_equal( sum( G[(nb1 + nb2 + 1):total_dim, (nb1 + nb2 + 1):total_dim]  - diag(ns)), 0)
})


testthat::test_that("get_smoothing_param_hybrid returns correct block-diagonal structure", {
  suppressPackageStartupMessages(library(fda))
  basis1 <- create.bspline.basis(c(0, 1), nbasis = 5)
  basis2 <- create.bspline.basis(c(0, 1), nbasis = 4)

  fd1 <- fd(coef = matrix(1, 5, 3), basisobj = basis1)
  fd2 <- fd(coef = matrix(2, 4, 3), basisobj = basis2)

  Z <- matrix(1, nrow = 3, ncol = 2)

  obj <- predictor_hybrid(Z = Z, functional_list = list(fd1, fd2))
  lambda <- c(0.5, 2)

  lambda_mat <- get_smoothing_param_hybrid(obj, lambda)

  # Check overall dimension
  expected_dim <- sum(obj$n_basis_list) + obj$n_scalar
  testthat::expect_equal(dim(lambda_mat), c(expected_dim, expected_dim))

  # Check diagonal entries
  testthat::expect_equal(Matrix::diag(lambda_mat)[1:5], rep(0.5, 5))
  testthat::expect_equal(Matrix::diag(lambda_mat)[6:9], rep(2, 4))
  testthat::expect_equal(Matrix::diag(lambda_mat)[10:11], rep(0, 2))  # scalar part
})


testthat::test_that("get_penalty_hybrid computes correct dimensions and structure", {
  suppressPackageStartupMessages(library(fda))

  # Setup
  basis <- fda::create.bspline.basis(c(0, 1), nbasis = 5)
  n_sample <- 3
  n_scalar <- 2

  fd1 <- fda::fd(coef = matrix(1, 5, n_sample), basisobj = basis)
  fd2 <- fda::fd(coef = matrix(2, 5, n_sample), basisobj = basis)
  Z <- matrix(rnorm(n_sample * n_scalar), n_sample, n_scalar)
  obj <- predictor_hybrid(Z = Z, functional_list = list(fd1, fd2))

  # Run
  penalty <- get_penalty_hybrid(obj)

  # Expected dimensions
  total_dim <- sum(obj$n_basis_list) + obj$n_scalar
  testthat::expect_equal(dim(penalty), c(total_dim, total_dim))

  # Check block diagonal structure: bottom-right should be zero matrix
  scalar_block <- as.matrix(penalty[
    (total_dim - n_scalar + 1):total_dim,
    (total_dim - n_scalar + 1):total_dim
  ])
  testthat::expect_equal(scalar_block, matrix(0, n_scalar, n_scalar))

  # Check that the upper blocks are positive semidefinite
  eigenvalues <- eigen(as.matrix(penalty[1:(total_dim - n_scalar), 1:(total_dim - n_scalar)]))$values
  testthat::expect_true(all(eigenvalues >= -1e-8))  # Numerical tolerance
})


testthat::test_that("get_pre_corrcov handles constant and random coefficients", {
  suppressPackageStartupMessages(library(fda))

  make_fd <- function(coefs, basis) {
    fd(coef = coefs, basisobj = basis)
  }

  basis <- create.bspline.basis(c(0, 1), nbasis = 4)
  n_sample <- 5
  n_scalar <- 2
  n_functional <- 2

  # Constant coefficients
  fd_const1 <- make_fd(matrix(1, 4, n_sample), basis)
  fd_const2 <- make_fd(matrix(2, 4, n_sample), basis)
  Z_const <- matrix(1, n_sample, n_scalar)
  y_const <- matrix(1:n_sample, ncol = 1)
  obj_const <- predictor_hybrid(Z = Z_const, functional_list = list(fd_const1, fd_const2))

  V_const <- get_pre_corrcov(obj_const, y_const)
  testthat::expect_equal(dim(V_const), rep(sum(obj_const$n_basis_list) + obj_const$n_scalar, 2))
  testthat::expect_true(isSymmetric(V_const))
  testthat::expect_true(all(V_const >= 0))

  # Random coefficients
  set.seed(123)
  fd_rand1 <- make_fd(matrix(rnorm(4 * n_sample), 4, n_sample), basis)
  fd_rand2 <- make_fd(matrix(rnorm(4 * n_sample), 4, n_sample), basis)
  Z_rand <- matrix(rnorm(n_sample * n_scalar), n_sample, n_scalar)
  y_rand <- matrix(rnorm(n_sample), ncol = 1)
  obj_rand <- predictor_hybrid(Z = Z_rand, functional_list = list(fd_rand1, fd_rand2))

  V_rand <- get_pre_corrcov(obj_rand, y_rand)
  testthat::expect_equal(dim(V_rand), rep(sum(obj_rand$n_basis_list) + obj_rand$n_scalar, 2))
  testthat::expect_true(isSymmetric(V_rand))
  testthat::expect_true(all(eigen(V_rand, only.values = TRUE)$values >= -1e-8))
})


