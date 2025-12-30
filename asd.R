

```




























### inprod.predictor_hybrid
- The inner product of $f_1=(f_1^{(1)}, \ldots,f_1^{(K)})$ and $f_2=(f_2^{(1)}, \ldots,f_2^{(K)})$ in $\mathcal{F}$ is defined as
$$\langle f_1, f_2\rangle_\mathcal{F} =  \sum_{k=1}^K \langle f_1^{(k)}, f_2^{(k)}\rangle_{L^2} = \sum_{k=1}^K \int_{\mathcal{T}_k} f_1^{(k)}(t_k) f_2^{(k)}(t_k) dt_k,$$ with norm
$$\Vert f_1 \Vert_\mathcal{F} = \langle f_1,f_1 \rangle_\mathcal{F}^{1/2} = \{ \sum_{k=1}^K \int_{\mathcal{T}_k} f_1^{(k)}(t_k)^2 dt_k\}^{1/2}.$$
  - We define the inner product between any two hybrid objects, $\mathbf{h}_1 = (f_1, \mathbf{v}_1)$ and $\mathbf{h}_2 = (f_2, \mathbf{v}_2)$, as
  $$
    \langle \mathbf{h}_1, \mathbf{h}_2\rangle_{\mathcal{H}} = \langle f_1, f_2\rangle_\mathcal{F} +  \langle  \mathbf{v_1}, \mathbf{v_2} \rangle = \sum \limits_{k=1}^K \int_{\mathcal{T}_k} f_1^{(k)}(t_k) f_2^{(k)}(t_k) dt_k + \sum \limits_{r=1}^p v_{1r}v_{2r},
    $$
      with norm $\Vert \cdot \Vert_\mathcal{H} = \langle \cdot,\cdot \rangle_\mathcal{H}^{1/2}$.
      
      - The method `inprod.predictor_hybrid` computes the inner product between two `predictor_hybrid` objects. It supports broadcasting when one of the inputs has a single observation.
      
      
      - Usage
      
      ```r
      inprod.predictor_hybrid(xi_1, xi_2)
      ```
      
      - Arguments:
        
        - `xi_1`: A `predictor_hybrid` object.
      - `xi_2`: Another `predictor_hybrid` object. If omitted, defaults to xi_1 (computes self-inner product).
      
      - Value:  A numeric vector of inner products (or a scalar if both inputs are single observations).
      
      - Details:
        - Functional components are summed using inprod() from the fda package.
      - Scalar components are handled via matrix multiplication.
      - Broadcasting is supported: if either xi_1 or xi_2 has only one sample, its values are broadcast across all rows of the other.
      - Ensures compatibility in the number of functional and scalar predictors before computing the result.
      ```{r}
      #' Inner product between two predictor_hybrid objects (with broadcasting)
      #'
      #' Computes the inner product between two `predictor_hybrid` objects,
      #' including both functional and scalar components. Supports broadcasting
      #' when one of the inputs has a single observation.
      #'
      #' @param xi_1 A `predictor_hybrid` object.
      #' @param xi_2 Another `predictor_hybrid` object. If missing, defaults to `xi_1`.
      #'
      #' @return A numeric vector of inner products, or a scalar if both inputs contain a single observation.
      #' @export
      inprod.predictor_hybrid <- function(xi_1, xi_2 = NULL) {
        # Handle self-inner product
        if (is.null(xi_2)) xi_2 <- xi_1
        
        # Type checks
        if (!inherits(xi_1, "predictor_hybrid") || !inherits(xi_2, "predictor_hybrid")) {
          stop("Both inputs must be of class 'predictor_hybrid'.")
        }
        
        # Structural checks
        if (xi_1$n_functional != xi_2$n_functional) {
          stop("Mismatch in number of functional predictors.")
        }
        if (xi_1$n_scalar != xi_2$n_scalar) {
          stop("Mismatch in number of scalar predictors.")
        }
        
        # Swap so that broadcasting always applies to xi_2
        if (xi_1$n_sample == 1 && xi_2$n_sample > 1) {
          tmp <- xi_1
          xi_1 <- xi_2
          xi_2 <- tmp
        }
        
        n1 <- xi_1$n_sample
        n2 <- xi_2$n_sample
        
        if (!(n1 == n2 || n2 == 1)) {
          stop("Sample sizes are incompatible for broadcasting.")
        }
        
        # Prepare components
        f1 <- xi_1$functional_list
        f2 <- xi_2$functional_list
        Z1 <- xi_1$Z
        Z2 <- xi_2$Z
        
        # Replicate fd and Z if needed
        if (n2 == 1) {
          f2 <- rep_fd(f2, n1)
          Z2 <- matrix(rep(c(Z2), n1), nrow = n1, byrow = TRUE)
        }
        
        # Compute functional inner products
        inprod_functional <- vapply(seq_len(n1), function(i) {
          sum(vapply(seq_along(f1), function(j) {
            fda::inprod(f1[[j]][i], f2[[j]][i])
          }, numeric(1)))
        }, numeric(1))
        
        
        # Compute scalar inner products
        inprod_scalar <- rowSums(Z1 * Z2)
        
        # Combine results
        result <- inprod_functional + inprod_scalar
        if (n1 == 1 && n2 == 1) as.numeric(result) else result
      }
      ```
      
      ```{r}
      #' Inner product between two predictor_hybrid objects (with broadcasting)
      #'
      #' Computes the inner product between two `predictor_hybrid` objects,
      #' including both functional and scalar components. Supports broadcasting
      #' when one of the inputs has a single observation.
      #'
      #' @param xi_1 A `predictor_hybrid` object.
      #' @param xi_2 Another `predictor_hybrid` object. If missing, defaults to `xi_1`.
      #'
      #' @return A numeric vector of inner products, or a scalar if both inputs contain a single observation.
      #' @export
      inprod_pen.predictor_hybrid <- function(xi_1, xi_2 = NULL, lambda) {
        
        inprod <- inprod.predictor_hybrid(xi_1, xi_2)
        # Prepare components
        
        f1 <- xi_1$functional_list
        
        f2 <- xi_2$functional_list
        for (j in 1:xi_1$n_functional){
          inprod <- inprod + lambda[j] * fda::inprod(
            fdobj1 = f1[[j]],
            fdobj2 = f2[[j]],
            Lfdobj1 = 2,
            Lfdobj2 = 2
          )
          
        }
        
        
        
        
        # Combine results
        return(as.numeric(inprod))
      }
      ```
      
      ### subset_predictor_hybrid
      
      ```{r}
      #' Extract a single observation from a predictor_hybrid object
      #'
      #' @param W A predictor_hybrid object with multiple samples.
      #' @param i Integer index of the sample to extract.
      #' @return A single-sample predictor_hybrid object.
      subset_predictor_hybrid <- function(W, i) {
        new_Z <- matrix(W$Z[i, ], nrow = 1)
        new_functional_list <- lapply(W$functional_list, function(fdobj) {
          fd(coef = matrix(coef(fdobj)[, i], ncol = 1), basisobj = fdobj$basis)
        })
        new_predictor <- predictor_hybrid(Z = new_Z, functional_list = new_functional_list)
        return(new_predictor)
      }
      ```
      
      ### replace_obs_hybrid
      
      ```{r}
      #' Replace a single observation in a predictor_hybrid object
      #'
      #' Replaces the i-th observation of a predictor_hybrid object with a new
      #' single-sample predictor_hybrid object.
      #'
      #' @param W A predictor_hybrid object with one or more samples.
      #' @param i An integer index specifying the observation to replace.
      #' @param new_W A single-sample predictor_hybrid object to use for replacement.
      #'
      #' @return A predictor_hybrid object with the i-th observation replaced.
      #' @export
      replace_obs_hybrid <- function(W, i, new_W) {
        # Input validation
        if (!inherits(W, "predictor_hybrid") || !inherits(new_W, "predictor_hybrid")) {
          stop("Both W and new_W must be of class 'predictor_hybrid'.")
        }
        if (new_W$n_sample != 1) {
          stop("new_W must be a single-sample predictor_hybrid object.")
        }
        if (i < 1 || i > W$n_sample) {
          stop(paste("Index i must be between 1 and", W$n_sample))
        }
        if (W$n_scalar != new_W$n_scalar) {
          stop("Mismatch in number of scalar predictors.")
        }
        if (W$n_functional != new_W$n_functional) {
          stop("Mismatch in number of functional predictors.")
        }
        
        # Check for compatible basis objects
        is_eqbasis <- getFromNamespace("is.eqbasis", "fda")
        for (j in seq_len(W$n_functional)) {
          if (!is_eqbasis(W$functional_list[[j]]$basis, new_W$functional_list[[j]]$basis)) {
            stop("Functional predictors must have the same basis objects.")
          }
        }
        
        # Replace the i-th observation in the scalar matrix Z
        W$Z[i, ] <- new_W$Z[1, ]
        
        # Replace the i-th observation in each functional predictor
        for (j in seq_along(W$functional_list)) {
          # The coefficients are stored as a matrix, with columns corresponding to samples
          W$functional_list[[j]]$coefs[, i] <- new_W$functional_list[[j]]$coefs[, 1]
        }
        
        # The number of samples remains the same
        return(W)
      }
      ```
      
      
      
      ### add_broadcast(for matrix)
      
      ```{r}
      #' @export
      #'
      #'
      # --- 1. Define the S4 Generic Function ---
      # This must be run first to tell R that 'add_broadcast' is a function
      # that can have methods defined for different classes.
      setGeneric("add_broadcast", function(input, other, alpha = 1) {
        standardGeneric("add_broadcast")
      })
      
      # --- 2. Define the S4 Method for 'matrix' Class ---
      setMethod("add_broadcast", "matrix",
                ###########################################
                function(input, other, alpha = 1){
                  
                  # 1. Structural Check: Ensure RHS is a single observation.
                  # This uses length() for vectors, and dim()[1] for matrices.
                  is_matrix_other <- is.matrix(other)
                  if (is_matrix_other && dim(other)[1] > 1) {
                    stop("RHS must have a single observation (row).")
                  }
                  
                  # 2. Conversion: Convert a single-row matrix into a vector for broadcasting.
                  # If 'other' is already a vector, this does nothing.
                  # If 'other' is a 1xP matrix, this flattens it into a length-P vector.
                  if (is_matrix_other) {
                    other <- as.vector(other)
                  }
                  
                  # 3. Native Arithmetic: Use simple R arithmetic for high performance.
                  # R broadcasts the vector 'other' across the rows of the matrix 'input'.
                  return( t( t(input) + alpha * other))
                }
                ###########################################
      )
      
      ```
      ### subtr_broadcast (for matrix)
      
      ```{r}
      #' @export
      #'
      #'
      setGeneric("subtr_broadcast", function(input, other, alpha = 1) {
        standardGeneric("subtr_broadcast")
      })
      setMethod("subtr_broadcast", "matrix",
                ###########################################
                function(input, other, alpha = 1){
                  return(add_broadcast(input, other, (-1*alpha)))
                }
                ###########################################
      )
      
      ```
      
      
      
      
      # One iteration
      
      ## small functions
      
      
      
      ## penalty matrix construction
      In this section, we provide an R function `get_constraint_matrix` that constructs the penalty matrix
      $$J^*+\Lambda \ddot{J}^\ast,$$
        as defined in Proposition 2 and used in the PLS component computation step of our proposed algorithm. We present functions  for computing $J^\ast$, $\Lambda$ and $\ddot{J}^\ast$ in sequence, and use these functions to define `get_constraint_matrix`.
      
      ### get_gram_matrix_block
      Constructs a block-diagonal Gram matrix for a hybrid predictor object, defined as
      $$
        J^*=\mathrm{blkdiag}(J, I_p) \in \mathbb{R}^{(MK+p) \times (MK+p)}
      $$
        where
      $J = \mathrm{blkdiag}(J^{(1)}, \cdots, J^{(K)}) \in \mathbb{R}^{MK \times MK}$, and
      $$
        J^{(k)} = \left[ \int_{\mathcal{T}k} b_m^{(k)}(t) , b_n^{(k)}(t) , dt \right]_{m,n=1}^M,
      $$
        as defined in \eqref{def:gram_basis} and \eqref{def:J_and_J_star}.
      
      - Arguments
      - obj: A `predictor_hybrid` object containing both functional and scalar components.
      
      - Value
      - A Matrix::bdiag sparse matrix representing the block-diagonal structure of the combined Gram matrix.
      
      - Details: This function builds a block-diagonal matrix by:
        - Stacking the Gram matrices of each functional component (from obj$gram_list),
      - Appending an identity matrix of size equal to the number of scalar predictors (to represent unpenalized scalar covariates).
      The resulting matrix has size $(\texttt{total_dim} \times \texttt{total_dim})$, where total_dim = \sum_k M_k + p, with $M_k$ the number of basis functions for the $k$-th functional predictor and $p$ the number of scalar covariates.
      
      - Usage
      ```r
      block_gram <- get_gram_matrix_block(my_predictor)
      ```
      
      
      **Code**
        ```{r}
      #' Construct block-diagonal Gram matrix for hybrid predictor
      #'
      #' Returns a block-diagonal matrix containing the Gram  matrices for
      #' each functional component and an identity matrix for the scalar part.
      #'
      #' @param obj A `predictor_hybrid` object.
      #'
      #' @return A block-diagonal matrix of size `(total_dim × total_dim)` where
      #' functional and scalar components are arranged in order.
      #' @export
      get_gram_matrix_block <- function(obj) {
        if (!inherits(obj, "predictor_hybrid")) {
          stop("Input must be of class 'predictor_hybrid'.")
        }
        
        gram_blocks <- c(obj$gram_list, list(diag(obj$n_scalar)))
        Matrix::bdiag(gram_blocks)
      }
      
      ```
      
      **Unit test**
        ```{r}
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
      
      ```
      ### get_smoothing_param_hybrid
      Constructs a block-diagonal smoothing parameter matrix for use in penalized estimation involving hybrid predictors, denoted $\Lambda \in \mathbb{R}^{(MK+p) \times (MK+p)}$,  defined in \eqref{def:Lambda} as
      $$
        \Lambda = \mathrm{blkdiag}(\lambda_1 I_M, \cdots, \lambda_K I_M, 0_{p \times p}) ,
      $$
        
        - Arguments
      - W: A predictor_hybrid object.
      - lambda: A numeric vector of length equal to the number of functional predictors. Each entry corresponds to a smoothing penalty weight.
      - Value: A sparse block-diagonal matrix combining scaled identity matrices for the functional parts and a zero matrix for the scalar part.
      - Details: This function generates a matrix used in regularized regression for functional predictors. Each functional component is penalized using a scaled identity matrix of its basis dimension. Scalar predictors are unpenalized.
      - Usage:
        ```r
      lambda_mat <- get_smoothing_param_hybrid(my_predictor, c(0.1, 0.2))
      ```
      **code**
        ```{r}
      #' Construct block-diagonal smoothing parameter matrix
      #'
      #' Generates a block-diagonal matrix with smoothing parameters applied to each
      #' functional component. Each block is a scaled identity matrix, where the scaling
      #' factor corresponds to the regularization parameter for that functional component.
      #' The scalar components are not penalized and thus contribute a zero matrix block.
      #'
      #' @param W A `predictor_hybrid` object.
      #' @param lambda A numeric vector of length equal to the number of functional components (`W$n_functional`), containing the smoothing parameters for each functional predictor.
      #'
      #' @return A block-diagonal matrix of size `(total_dim × total_dim)`, where the top-left blocks are scaled identity matrices for functional predictors and the bottom-right block is a zero matrix for scalar covariates.
      #' @export
      get_smoothing_param_hybrid <- function(W, lambda) {
        if (!inherits(W, "predictor_hybrid")) {
          stop("Input W must be of class 'predictor_hybrid'.")
        }
        
        if (length(lambda) != W$n_functional) {
          stop("Length of lambda must match the number of functional predictors.")
        }
        
        lambda_blocks <- lapply(seq_len(W$n_functional), function(ii) {
          nb <- W$functional_list[[ii]]$basis$nbasis
          lambda[ii] * diag(nb)
        })
        
        lambda_blocks[[W$n_functional + 1]] <- matrix(0, nrow = W$n_scalar, ncol = W$n_scalar)
        
        Matrix::bdiag(lambda_blocks)
      }
      ```
      
      **unit test**
        ```{r}
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
      
      ```
      
      
      ### get_penalty_hybrid
      Constructs a matrix defined in \label{def:J_dotdot_ast}
      $$
        \ddot{J}^\ast=\mathrm{blkdiag}(\ddot{J}^{(1)}, \cdots, \ddot{J}^{(K)}, 0_{p \times p}) \in \mathbb{R}^{(MK+p) \times (MK+p)},
      $$
        where
      $\ddot{J}^{(k)}$ is  a $M \times M$ matrix defined in \label{def:dodot_J_k} as the gram matrix formed by the second derivativ of the basis functions
      $$
        \ddot{J}^{(k)} = \left[ \int_{\mathcal{T}_k} \ddot{b}^{(k)}_m(t) \ddot{b}^{(k)}_n(t) \, dt \right]_{m,n=1}^M.
      $$
        
        
        
        ```{r}
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
      
      ```
      
      ### get_constraint_matrix
      
      ```r
      W <- predictor_hybrid(Z = Z, functional_list = list(fd1, fd2))
      lambda <- c(0.1, 0.2)
      constraint_matrix <- get_constraint_matrix_hybrid(W, lambda)
      ```
      
      **code**
        ```{r}
      #' Construct penalized constraint matrix for hybrid predictors
      #'
      #' Computes the denominator matrix used in penalized estimation for hybrid predictors,
      #' combining the Gram matrix, smoothing parameter matrix, and penalty matrix.
      #'
      #' Specifically, this function returns the matrix:
      #' \deqn{J^\ast + \Lambda \ddot{J}^\ast}
      #' where \eqn{J^\ast} is the block-diagonal Gram matrix,
      #' \eqn{\Lambda} is the block-diagonal smoothing parameter matrix,
      #' and \eqn{\ddot{J}^\ast} is the block-diagonal penalty matrix of second derivative inner products.
      #'
      #' @param W A `predictor_hybrid` object.
      #' @param lambda A numeric vector of smoothing parameters, one for each functional predictor.
      #'
      #' @return A matrix representing the penalized constraint matrix used in estimation.
      #' @export
      get_constraint_matrix_hybrid <- function(W, lambda) {
        J_star <- get_gram_matrix_block(W)
        J_dotdot_star <- get_penalty_hybrid(W)
        lambda_mat <- get_smoothing_param_hybrid(W, lambda)
        
        J_star + lambda_mat %*% J_dotdot_star
      }
      
      ```
      
      ## PLS component computation
      
      ### Spectral method
      The old implementation that computes the leading eigenvector of
      
      ### helper: get_CJ
      ```{r}
      get_CJ <- function(W){
        n <- W$n_sample
        K <- W$n_functional
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
        return(C_J)
      }
      ```
      
      
      
      ### helper: get_pre_corrcov
      construct the matrix V_star:
        $$
        V^\ast = n^{-2}\begin{bmatrix}
      J \widetilde{C}^\top \mathbf{y}\mathbf{y}^\top \widetilde{C} J & J \widetilde{C}^\top \widetilde{\mathbf{y}} \widetilde{\mathbf{y}}^\top \widetilde{Z} \\
      \widetilde{Z}^\top \widetilde{\mathbf{y}} \widetilde{\mathbf{y}}^\top \widetilde{C} J &   \widetilde{Z}^\top \widetilde{\mathbf{y}} \widetilde{\mathbf{y}}^\top \widetilde{Z}
      \end{bmatrix} = n^{-2} J^* \widetilde{C}^{*\top} \widetilde{\mathbf{y}} \widetilde{\mathbf{y}}^\top \widetilde{C}^{*} J^*
        \in \mathbb{R}^{(MK+p) \times (MK+p)},
      $$
        where $J$ is defined in `get_gram_matrix_block`
      equation (7),
      ```{r}
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
        C_J <- get_CJ(W)
        
        
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
      ```
      
      
      
      
      ### main function (wrapper: get_pls_comp
      Computes the first PLS component from a hybrid predictor and response vector using regularization.
      
      - Inputs
      - **`W`**: A `predictor_hybrid` object with functional and scalar predictors.
      - **`y`**: A numeric response vector of length equal to `W$n_sample`.
      - **`L`**: A Cholesky decomposition of a positive definite regularization matrix.
      - Output:
        - `xi_hat`: The first PLS component as a `predictor_hybrid` object
      - `E`: The matrix `inv(L) %*% V_star %*% t(inv(L))`
      - `V_star`: Cross-covariance matrix between predictors and response
      - `eigen_val`: Leading eigenvalue of `E`
      ```{r}
      #' Compute the First PLS Component from a Hybrid Predictor
      #'
      #' Computes the coefficients of the first Partial Least Squares (PLS) component based on
      #' a hybrid predictor object and a response vector, using a regularized generalized eigenvalue problem.
      #'
      #' @param W A predictor_hybrid object containing both functional and scalar predictors.
      #' @param y A numeric response vector of length equal to the number of samples in W.
      #' @param L cholesky decompisition of a regularization matrix (typically positive definite) used in the generalized eigenproblem.
      #'
      #' @return A list with the following elements:
      #'
      #'   xi_hat: The estimated first PLS component as a predictor_hybrid object.
      #'   E: The generalized eigenvalue problem matrix E = inv(L) V* t(inv(L)).
      #'   V_star: The cross-covariance matrix between predictors and response.}
      #'   eigen_val: The leading eigenvalue of E.
      #'
      #'
      #' @export
      get_pls_comp <- function(W, y, L){
        
        
        V_star <- get_pre_corrcov(W, y)
        
        #invL <- Matrix::chol2inv(L)
        invL <- solve(L)
        A <- V_star %*% t(invL)       # A = V* × t(inv(L))
        E <- invL %*% A               # E = inv(L) × A
        
        eigen_result <- eigen(E)
        e <- eigen_result$vectors[, 1]
        if (is.complex(e)){
          print("stop")
          return("stop")
        }
        
        
        xi_star <- t(invL) %*% e      # xi_star solves t(L) xi = e
        xi_hat <- predictor_hybrid_from_coef(format = W, coef = xi_star)
        return(xi_hat)
      }
      ```
      
      ## linear system method
      
      ### helper: get_xi_hat_linear
      ```{r}
      get_xi_hat_linear_pen <- function(W, y, lambda){
        # new
        ## step 1: inner products
        n <- W$n_sample
        K <- W$n_functional
        u <- gamma <- list()
        
        v <- t(W$Z) %*% y
        q <- sum(v^2)
        for (j in 1:K){
          Theta_t <- W$functional_list[[j]]$coefs
          B <- W$gram_list[[j]]
          u[[j]] <- B %*% Theta_t %*% y
          
          
          # linear system
          R <- W$gram_list[[j]] + lambda[j] * W$gram_deriv2_list[[j]]
          gamma[[j]] <- solve(R, u[[j]])
          
          #q
          q <- q + sum( gamma[[j]] * ( R %*% gamma[[j]] ))
        }
        d_vec <- c(do.call(c, gamma), v)  # or: unlist(d)
        ## step 2: squared norm
        d_vec <- d_vec/ sqrt(q)
        xi_hat <- predictor_hybrid_from_coef(format = W, coef = d_vec)
        return(xi_hat)
      }
      ```
      
      
      the full vector of PLS scores $$\hat{\boldsymbol{\rho}}  = (\hat{\rho}_1, \ldots, \hat{\rho}_n)^\top$$
        is computed through the following matrix multiplication:
        $$
        \hat{\boldsymbol{\rho}}^\top = \sum_{k=1}^K ( \hat{\boldsymbol{\gamma}}_j )^\top B \, \Theta_j^\top + \boldsymbol{\zeta}^\top \mathbf{Z}^\top.
      $$
        ```{r}
      get_rho <- function(d_vec, W){
        K <- W$n_functional
        M <- W$n_basis_list[[1]]
        n_scalar <- W$n_scalar
        n <- W$n_sample
        rho <- rep(0, n)
        
        #functional part
        for (j in 1:K){
          start_index <- (j-1)*M + 1
          end_index <- j*M
          gamma_j <- d_vec[ start_index : end_index]
          Theta_j  <- t(W$functional_list[[j]]$coefs)
          B <- W$gram_list[[j]]
          rho <- rho + Theta_j %*% B %*% gamma_j
        }
        
        #scalar part
        rho <- rho + W$Z %*% d_vec[(K * M + 1):length(d_vec)]
        return (as.vector(rho))
      }
      ```
      
      
      ```{r}
      library(fda)
      
      # --- Step 1: Simulate Data ---
      set.seed(111)
      
      n_sample <- 100
      n_scalar <- 2
      n_basis <- 20
      n_functional <- 2
      
      # Scalar predictors Z (centered)
      Z <- matrix(rnorm(n_sample * n_scalar), ncol = n_scalar)
      Z <- scale(Z, center = TRUE, scale = FALSE)
      
      # Create functional predictors (centered)
      basis <- create.bspline.basis(c(0, 1), nbasis = n_basis)
      functional_list <- list()
      for (j in 1:n_functional) {
        coefs <- matrix(rnorm(n_basis * n_sample), nrow = n_basis)
        fdobj <- fd(coef = coefs, basisobj = basis)
        mean_fd <- mean.fd(fdobj)
        mean_coef <- coef(mean_fd)
        centered_coefs <- coefs - matrix(mean_coef, nrow = n_basis, ncol = n_sample)
        functional_list[[j]] <- fd(coef = centered_coefs, basisobj = basis)
      }
      
      # --- Step 2: Construct hybrid object ---
      W <- predictor_hybrid(Z = Z, functional_list = functional_list)
      
      # --- Step 3: Simulate and center response ---
      total_coef_len <- sum(W$n_basis_list) + W$n_scalar
      true_coef <- rnorm(total_coef_len)
      
      # Flatten design matrix
      Xmat <- matrix(NA, n_sample, total_coef_len)
      for (i in 1:n_sample) {
        W_i <- subset_predictor_hybrid(W, i)
        func_coefs <- do.call(c, lapply(W_i$functional_list, function(fdobj) as.vector(coef(fdobj))))
        Xmat[i, ] <- c(func_coefs, W_i$Z)
      }
      
      y <- Xmat %*% true_coef + rnorm(n_sample, sd = 0.1)
      y <- as.vector(scale(y, center = TRUE, scale = FALSE))
      
      # --- Step 4: Compute constraint matrix ---
      lambda <- rep(0, W$n_functional)
      #constr_mat <- get_constraint_matrix_hybrid(W, lambda)
      #L <- chol(constr_mat)
      
      # --- Step 5: Run eigen-based PLS ---
      W_template <- W  # prevent mutation
      
      
      pls_eig_result <- get_pls_comp(W, y, L)
      xi_hat_eig <- pls_eig_result$xi_hat
      
      xi_hat_pen <- get_xi_hat_linear_pen(W, y, lambda)
      xi_hat_lin_new <- get_xi_hat_linear_new(W, y)
      # Wrap into predictor_hybrid object using format template
      xi_hat_direct <- predictor_hybrid_from_coef(format = W_template, coef = xi_hat_lin)
      norm_val <- sqrt(inprod.predictor_hybrid(xi_hat_direct))
      xi_hat_direct <- scalar_mul.predictor_hybrid(xi_hat_direct,  1 / norm_val)
      xi_hat_direct_new <- predictor_hybrid_from_coef(format = W_template, coef = xi_hat_lin_new)
      
      cat(xi_hat_direct_new$functional_list[[1]]$coefs- xi_hat_direct$functional_list[[1]]$coefs)
      cat(xi_hat_direct_new$functional_list[[2]]$coefs- xi_hat_direct$functional_list[[2]]$coefs)
      
      
      ```
      
      
      
      ## Residualization
      
      ### response residualization
      
      ```{r}
      get_nu <- function(y, rho){
        nu <- sum(y*rho) / sum(rho*rho)
        return(nu)
      }
      ```
      
      ```{r}
      residualize_y <- function(y, rho, nu){
        y_next <- y - nu * rho
        return(y_next)
      }
      ```
      
      
      
      
      ```{r}
      LSE_ptws <- function(input, rho) {
        # old implementation
        rho_norm_sq <- (Matrix::norm(rho, type = "2"))^2
        nu <- as.numeric(t(rho) %*% input) / rho_norm_sq
        return(nu)
      }
      
      ```
      
      ###  predictor residualization
      
      $$
        \widetilde{W}_i^{[l+1]}
      := 	\widetilde{W}_i^{[l]}  -
        \widehat{\rho}_i^{[l]}
      \hat{\delta}^{[l]},~\text{where}~
        \delta^{[l]}
      :=
        \frac{1}{\| 	\hat{\boldsymbol{\rho}}^{[l]}\|_2^2}
      \sum_{i=1}^n
      \widehat{\rho}_i^{[l]}
      \widetilde{W}_i^{[l]},
      $$
        
        
        ```{r}
      get_delta <- function(W, rho){
        delta <- scalar_mul.predictor_hybrid(subset_predictor_hybrid(W, 1), rho[1])
        for (i in 2:length(rho)) {
          delta <- add.predictor_hybrid(delta, subset_predictor_hybrid(W, i),  rho[i])
        }
        delta <- scalar_mul.predictor_hybrid(delta, 1/sum(rho*rho))
        return(delta)
      }
      ```
      
      
      ```{r}
      residualize_predictor <- function(W, rho, delta){
        n <- length(rho)
        W_res <- W
        for (i in 1:n){
          W_i <- subset_predictor_hybrid(W, i)
          W_res_i <-subtr.predictor_hybrid(W_i, delta,  rho[i])
          W_res <- replace_obs_hybrid(W_res, i, W_res_i)
        }
        return(W_res)
      }
      ```
      
      
      
      
      ### LSE_hybrid
      
      ```{r}
      LSE_hybrid <- function(W, rho) {
        # old implementation
        C_star_J_star <- matrix(NA, nrow = W$n_sample, ncol = sum(W$n_basis_list) + W$n_scalar)
        for (ii in (1 : W$n_functional)){
          cumsum_n_basis_now <- sum(W$n_basis_list[1:ii])
          n_basis_now <- W$n_basis_list[ii]
          C_star_J_star[, (cumsum_n_basis_now - n_basis_now + 1) :cumsum_n_basis_now] <-
            t(W$functional_list[[ii]]$coefs) %*% W$jacobian_list[[ii]]
        }
        C_star_J_star[ ,(sum(W$n_basis_list) + 1) : ncol(C_star_J_star)] <- W$Z
        rho_t_C_star_J_star <- t(rho) %*% C_star_J_star
        
        # matrices in the statement of Proposition 4
        J_star <-as.matrix(get_jacobian_hybrid(W))
        d_hat_star <- t( #basis coefficient of hybrid regression coefficientdelta
          solve(
            t((norm(rho, type="2"))^2 * J_star),
            t(rho_t_C_star_J_star)
          )
        )
        delta_hat <- hybrid_from_coef(W, d_hat_star) #hybrid regression coefficient
        return(delta_hat)
      }
      ```
      
      
      
      # main algorithm
      
      ```{r eval=FALSE, include=FALSE}
      #' @export
      fit.hybridPLS <- function(W, y, n_iter, lambda) {
        # 1. Initialize storage
        W_now <- rho <- xi <- delta <- nu  <- iota <- beta <- list()
        
        # 3. Initialize residuals
        W_now[[1]] <- W
        y_now <- y
        for (l in 1:n_iter) {
          cat(paste(l, "th component", "\n"))
          
          xi[[l]]  <- get_xi_hat_linear_pen(W_now[[l]], y_now, lambda) # PLS direction
          rho[[l]] <- inprod.predictor_hybrid(W_now[[l]], xi[[l]]) #PLS score
          delta[[l]] <- get_delta(W_now[[l]], rho[[l]]);
          W_now[[l + 1]] <- residualize_predictor(W_now[[l]], rho[[l]], delta[[l]]) # predictor residual
          nu[[l]] <- get_nu(y_now, rho[[l]]); y_now <- residualize_y(y_now, rho[[l]], nu[[l]]) # response residual
          
          iota[[l]] <- xi[[l]]
          if (l == 1) {
            beta[[l]] <- scalar_mul.predictor_hybrid( iota[[l]], nu[[l]] )
          } else {
            for (u in 1:(l - 1)) iota[[l]] <- subtr.predictor_hybrid( iota[[l]], iota[[u]], inprod.predictor_hybrid(delta[[u]], xi[[l]]))
            beta[[l]] <- add.predictor_hybrid(beta[[l - 1]], iota[[l]], nu[[l]])
          }
        }
        return (list(
          rho = rho,
          xi = xi,
          W = W_now,
          beta = beta
        ))
      }
      
      ```
      
      
      ```{r eval=FALSE, include=FALSE}
      #' @export
      fit.hybridPLS_eigen <- function(W, y, n_iter, lambda) {
        # 1. Initialize storage
        W_now <- rho <- xi <- delta <- nu  <- sigma <- eta <- list()
        
        # 3. Initialize residuals
        W_now[[1]] <- W
        y_now <- y
        final_succesful_iteration <- 0
        constr_mat <- get_constraint_matrix_hybrid(W, lambda)
        constr_mat_chol <- t(chol(constr_mat))
        for (l in 1:n_iter) {
          cat(paste(l, "th component", "\n"))
          
          #xi[[l]]  <- get_xi_hat_linear_pen(W_now[[l]], y_now, lambda) # PLS direction
          xi[[l]]  <- get_pls_comp(W_now[[l]], y_now, constr_mat_chol) # PLS direction
          rho[[l]] <- inprod.predictor_hybrid(W_now[[l]], xi[[l]]) #PLS score
          delta[[l]] <- get_delta(W_now[[l]], rho[[l]]);
          W_now[[l + 1]] <- residualize_predictor(W_now[[l]], rho[[l]], delta[[l]]) # predictor residual
          nu[[l]] <- get_nu(y_now, rho[[l]]); y_now <- residualize_y(y_now, rho[[l]], nu[[l]]) # response residual
          
          # Step 3: Orthogonalize and update prediction
          #sigma[[l]] <- xi[[l]]
          #if (l == 1) {
          #  eta[[l]] <- scalar_mul(sigma[[l]], nu[[l]])
          #} else {
          #  for (u in 1:(l - 1)) {
          #    sigma[[l]] <- subtr(sigma[[l]], sigma[[u]], hybrid_inner_prod(delta[[u]], xi[[l]]))
          #  }
          #  eta[[l]] <- add(eta[[l - 1]], sigma[[l]], nu[[l]])
          #}
        }
        return (list(
          rho = rho,
          xi = xi,
          W = W_now
        ))
      }
      
      ```
      
      ```{r}
      inprod.predictor_hybrid(res$W[[5]], res$xi[[3]])
      
      sum(y * inprod.predictor_hybrid(res$W[[5]], res$xi[[3]]))
      
      inprod_pen.predictor_hybrid(res$xi[[1]], res$xi[[8]], lambda) - inprod_pen.predictor_hybrid(res_eigen$xi[[1]], res_eigen$xi[[8]], lambda)
      inprod_pen.predictor_hybrid(res$xi[[1]], res$xi[[8]], 10*lambda) - inprod_pen.predictor_hybrid(res_eigen$xi[[1]], res_eigen$xi[[8]], 10*lambda)
      inprod_pen.predictor_hybrid(res$xi[[1]], res$xi[[8]], lambda/10) - inprod_pen.predictor_hybrid(res_eigen$xi[[1]], res_eigen$xi[[8]], lambda/10)
      sum(y * inprod.predictor_hybrid(res$W[[8]], res$xi[[1]])) - sum(y * inprod.predictor_hybrid(res_eigen$W[[8]], res_eigen$xi[[1]]))
      
      
      
      
      
      
      inprod_pen.predictor_hybrid(res$xi[[3]], res$xi[[4]], lambda) - inprod_pen.predictor_hybrid(res_eigen$xi[[3]], res_eigen$xi[[4]], lambda)
      inprod_pen.predictor_hybrid(res$xi[[3]], res$xi[[4]], 10*lambda) - inprod_pen.predictor_hybrid(res_eigen$xi[[3]], res_eigen$xi[[4]], 10*lambda)
      inprod_pen.predictor_hybrid(res$xi[[3]], res$xi[[4]], lambda/10)
      sum(y * inprod.predictor_hybrid(res$W[[3]], res$xi[[1]]))
      
      
      inprod_pen.predictor_hybrid(res$xi[[5]], res$xi[[10]], lambda)
      inprod_pen.predictor_hybrid(res$xi[[5]], res$xi[[10]], 10*lambda)
      inprod_pen.predictor_hybrid(res$xi[[5]], res$xi[[10]], lambda/10)
      sum(y * inprod.predictor_hybrid(res$W[[5]], res$xi[[1]]))
      
      inprod_pen.predictor_hybrid(res$xi[[4]], res$xi[[4]], lambda)
      inprod_pen.predictor_hybrid(res$xi[[4]], res$xi[[4]], lambda/10)
      inprod_pen.predictor_hybrid(res$xi[[4]], res$xi[[4]], 10*lambda )
      
      ```
      
      
      
      
      ```{r}
      lambda <- c(0.5,0.5)
      L=8
      res <- fit.hybridPLS(W,y,L, lambda)
      res_eigen <- fit.hybridPLS_eigen(W,y,L, lambda)
      ```
      
      
      ### nipals_pen_hybrid
      
      
      - depends on:
        - 'get_constraint_matrix_hybrid'
      - 'get_pls_comp'
      - 'hybrid_inner_prod'
      - 'LSE_hybrid'
      
      ```{r eval=FALSE, include=FALSE}
      #' Iterative Penalized NIPALS Algorithm for Hybrid PLS
      #'
      #' Performs iterative extraction of PLS components using a penalized NIPALS algorithm
      #' for hybrid predictor objects that include both scalar and functional data.
      #'
      #' At each iteration:
      #' \enumerate{
      #'   \item Computes the first penalized PLS component using \code{get_pls_comp}.
      #'   \item Computes PLS scores via hybrid inner products.
      #'   \item Residualizes both predictors and response with respect to the extracted scores.
      #'   \item Orthogonalizes and accumulates loading vectors to construct final linear predictors.
      #' }
      #' The same constraint matrix (computed once) is used across all iterations.
      #'
      #' @param W A \code{predictor_hybrid} object containing scalar and functional predictors.
      #' @param y A numeric response vector of length equal to \code{W$n_sample}.
      #' @param n_iter An integer specifying the number of PLS components (iterations) to compute.
      #' @param lambda A numeric vector of smoothing parameters (length = number of functional predictors).
      #' @param verbose Logical; if \code{TRUE}, prints progress for each component.
      #'
      #' @return This function operates by side effects. Internally, it constructs and updates:
      #' \itemize{
      #'   \item \code{xi}, \code{rho}, \code{delta}, \code{nu}: component vectors and regression effects,
      #'   \item \code{sigma}, \code{eta}: orthogonalized components and accumulated prediction functions,
      #'   \item \code{W_now}, \code{resid_y}: residualized predictor and response at each iteration,
      #'   \item \code{E}, \code{V_star}, \code{eigen_val}: matrices for monitoring stability and eigenvalue progress.
      #' }
      #' The final number of successful iterations is tracked via \code{final_succesful_iteration}.
      #'
      #' @export
      nipals_pen_hybrid <- function(W, y, n_iter, lambda, verbose) {
        # 1. Initialize storage
        rho <- xi <- delta <- nu  <- sigma <- eta <- list()
        E <- V_star <- eigen_val <- list()
        fitted_value_W <- fitted_value_y <- list()
        resid_y <- W_now <- list()
        first_eigen_val <- mse_W <- mse_y <- rep(NA, n_iter)
        
        # 2. Compute and store the constraint matrix once
        constr_mat <- get_constraint_matrix_hybrid(W, lambda)
        constr_mat_chol <- t(chol(constr_mat))
        
        # 3. Initialize residuals
        W_now[[1]] <- W
        y_now <- y
        final_succesful_iteration <- 0
        
        for (l in 1:n_iter) {
          if (verbose) cat(paste(l, "th component", "\n"))
          
          # Step 1: Compute penalized PLS direction and score
          pls_result_now <- get_pls_comp(W_now[[l]], y_now, constr_mat_chol, verbose)
          if (!is.list(pls_result_now)) break
          
          xi[[l]] <- pls_result_now$xi
          rho[[l]] <- hybrid_inner_prod(W_now[[l]], xi[[l]])
          
          # Step 2: Residualize predictor and response
          delta[[l]] <- LSE_hybrid(W_now[[l]], rho[[l]])
          fitted_value_W[[l]] <- fitted_value(delta[[l]], rho[[l]])
          W_now[[l + 1]] <- subtr(W_now[[l]], fitted_value_W[[l]])
          
          nu[[l]] <- LSE_ptws(y_now, rho[[l]])
          fitted_value_y[[l]] <- nu[[l]] * rho[[l]]
          y_now <- y_now - fitted_value_y[[l]]
          resid_y[[l]] <- y_now
          
          # Step 3: Orthogonalize and update prediction
          sigma[[l]] <- xi[[l]]
          if (l == 1) {
            eta[[l]] <- scalar_mul(sigma[[l]], nu[[l]])
          } else {
            for (u in 1:(l - 1)) {
              sigma[[l]] <- subtr(sigma[[l]], sigma[[u]], hybrid_inner_prod(delta[[u]], xi[[l]]))
            }
            eta[[l]] <- add(eta[[l - 1]], sigma[[l]], nu[[l]])
          }
          
          # Record diagnostics
          E[[l]] <- pls_result_now$E
          V_star[[l]] <- pls_result_now$V_star
          eigen_val[[l]] <- pls_result_now$eigen_val
          first_eigen_val[l] <- eigen_val[[l]][1]
          final_succesful_iteration <- l
        }
      }
      
      ```
      
      
      # Simulation tools
      
      ## sample splitting and normalization
      
      ```{r}
      create_idx_train_test <-function(n_sample, train_ratio){
        n_train <- floor(n_sample*train_ratio)
        idx_train <- sort(sample(n_sample,n_train, replace=F))
        idx_test <- sort((1 : n_sample)[-idx_train])
        return(list(
          idx_train = idx_train, idx_test = idx_test
        ))
      }
      ```
      
      
      ```{r}
      #' @title Generate Train Indices
      #' @description Generates train/test split indices.
      #' @param n_sample Integer. Number of samples.
      #' @param train_ratio Numeric. Ratio of train samples.
      #' @return Vector of train indices.
      get_idx_train <- function(n_sample, train_ratio) {
        n_train <- floor(n_sample * train_ratio)
        idx_train <- sort(sample(n_sample, n_train, replace = FALSE))
        return(idx_train)
      }
      ```
      
      
      ```{r}
      create_idx_kfold <- function(n_sample, n_fold){
        idx_list <- caret::createFolds( (1 : n_sample), k = n_fold, list = TRUE, returnTrain = FALSE)
        idx_list_list <- list()
        for (i in 1:length(idx_list)){
          valid_hyper_index <- idx_list[[i]]
          train_hyper_index <- (1 : n_sample)[-valid_hyper_index]
          idx_list_list[[i]] <- list("idx_train" = train_hyper_index, "idx_valid" = valid_hyper_index)
        }
        return(idx_list_list)
      }
      ```
      
      
      
      ```{r}
      n_sample.fd <- function(fd_obj){
        length(fd_obj$fdnames$reps)
      }
      ```
      
      
      ### split.all
      
      Splits functional predictors, scalar predictors, and response into train and test sets.
      
      ```{r}
      #' Split Hybrid Predictor and Response Data
      #'
      #' Splits the hybrid predictor object and the response vector into training and
      #' testing sets based on a given training ratio.
      #'
      #' @param W_hybrid predictor_hybrid object containing all functional and scalar data.
      #' @param response Vector. The response variable corresponding to the samples in W_hybrid.
      #' @param train_ratio Numeric scalar between 0 and 1 indicating the proportion of data to use for training.
      #' @return List containing split train and test sets for all components.
      #' @export
      split.all <- function(W_hybrid, response, train_ratio) {
        # --- Type and Ratio Checks ---
        stopifnot(
          "'W_hybrid' must be a predictor_hybrid object" = inherits(W_hybrid, "predictor_hybrid"),
          "'response' must be a vector" = is.vector(response),
          "Response length must match hybrid sample size" = length(response) == W_hybrid$n_sample,
          "'train_ratio' must be between 0 and 1" = is.numeric(train_ratio) && train_ratio > 0 && train_ratio < 1
        )
        
        N <- W_hybrid$n_sample
        N_train <- floor(N * train_ratio)
        
        # Generate random indices for training set
        # Using 1:N ensures indices are always valid
        train_idx <- sample(1:N, N_train)
        test_idx <- (1:N)[-train_idx]
        
        W_train <- predictor_hybrid(
          W_hybrid$Z[train_idx, , drop = FALSE],
          list(
            W_hybrid$functional_list[[1]][train_idx],
            W_hybrid$functional_list[[2]][train_idx]
          )
        )
        
        W_test <- predictor_hybrid(
          W_hybrid$Z[test_idx, , drop = FALSE],
          list(
            data_fun_1_test = W_hybrid$functional_list[[1]][test_idx],
            data_fun_2_test = W_hybrid$functional_list[[2]][test_idx]
          )
        )
        
        return(
          list(
            predictor_train = W_train,
            predictor_test = W_test ,
            response_train = response[train_idx],
            response_test = response[test_idx]
          ))
      }
      ```
      
      
      
      
      ### curve_normalize
      Subtracts the mean function and divides by a scaling constant (deno). This implementation does not use the native fda S3 methods (minus.fd and times.fd), since they do not broadcast well.
      ```{r}
      #'
      #' @export
      curve_normalize <- function(functional, mean_functional, deno){
        # operation also takes care of the original_X
        
        functional_normalized <- fd(coef =
                                      coef(functional) - coef(mean_functional) %*% matrix(1, ncol = ncol(coef(functional))),
                                    basisobj = functional$basis
        )#numerator
        functional_normalized <- times.fd(1/deno, functional_normalized)
        return(functional_normalized)
      }
      ```
      
      ### curve_normalize_train_test
      Normalizes the functional components of both the training and testing 'predictor_hybrid' objects using the mean and standard deviation derived ONLY from the training data. The scalar components (Z) are left unchanged.
      
      ```{r}
      #' Normalize Functional Predictors in Hybrid Objects (Train/Test Split)
      #' @param train A `predictor_hybrid` object representing the training set.
      #' @param test A `predictor_hybrid` object representing the testing set.
      #' @return A list containing:
      #'   train: The normalized training `predictor_hybrid` object.}
      #'   test: The normalized testing `predictor_hybrid` object.}
      #'   deno_train: A numeric vector of scaling factors used for each functional predictor.}
      #' @export
      curve_normalize_train_test <- function(train, test) {
        # --- Type Checks ---
        stopifnot(
          "'train' must be of class 'predictor_hybrid'" = inherits(train, "predictor_hybrid"),
          "'test' must be of class 'predictor_hybrid'" = inherits(test, "predictor_hybrid"),
          "Mismatch in number of functional predictors" = train$n_functional == test$n_functional
        )
        
        train_normalized <- train
        test_normalized <- test
        deno_train <- numeric(train$n_functional)
        mean_train <- list()
        
        # Iterate through each functional predictor (k)
        for (k in seq_len(train$n_functional)) {
          # 1. Determine centering and scaling factors from the TRAINING data only
          train_fd <- train$functional_list[[k]]
          
          # Calculate the mean function for centering
          mean_train[[k]] <- fda::mean.fd(train_fd)
          
          # Calculate the standard deviation function
          sd_train <- fda::sd.fd(train_fd)
          
          # Calculate the scaling factor (L2 norm of the standard deviation function)
          deno_train[k] <- sqrt(fda::inprod(sd_train, sd_train))
          
          # Ensure scaling factor is positive
          if (deno_train[k] <= 1e-10) {
            warning(paste("Scaling factor for functional predictor", k, "is near zero. Skipping normalization for this predictor."))
            deno_train[k] <- 1 # Set to 1 to prevent division by zero, effectively no normalization
          }
          
          # 2. Normalize the TRAINING set
          train_normalized$functional_list[[k]] <- curve_normalize(
            train_fd, # Use original train data
            mean_train[[k]],
            deno_train[k]
          )
          
          # 3. Normalize the TEST set (using train factors to prevent leakage)
          test_normalized$functional_list[[k]] <-
            curve_normalize(
              test$functional_list[[k]], # Use original test data
              mean_train[[k]], # Centered by training mean
              deno_train[k] # Scaled by training SD
            )
        }
        # Return the modified hybrid objects and the scaling factors
        return( list(
          predictor_train = train_normalized,
          predictor_test = test_normalized,
          mean_train = mean_train,
          deno_train = deno_train
        ))
      }
      ```
      
      
      ```{r}
      #'
      #' @export
      scalar_normalize <- function(scalar_predictor, mean, deno){
        scalar_predictor_normalized <- subtr_broadcast(scalar_predictor, mean) #numerator
        deno_mat <- matrix(rep(1,nrow(scalar_predictor))) %*% deno
        scalar_predictor_normalized <- scalar_predictor_normalized / deno_mat
        return(scalar_predictor_normalized)
      }
      ```
      
      
      ```{r}
      #' Normalize the scalar component (Z) of hybrid predictors using training set statistics.
      #'
      #' @param predictor_train A `predictor_hybrid` object (train set).
      #' @param predictor_test A `predictor_hybrid` object (test set).
      #' @return A list containing the normalized 'train' and 'test' hybrid objects.
      #' @export
      scalar_normalize_train_test <- function(predictor_train, predictor_test){
        
        # Type and structure checks
        stopifnot(
          inherits(predictor_train, "predictor_hybrid"),
          inherits(predictor_test, "predictor_hybrid"),
          predictor_train$n_scalar == predictor_test$n_scalar
        )
        
        # Create copies for modification
        predictor_train_normalized <- predictor_train
        predictor_test_normalized <- predictor_test
        
        # 1. Calculate statistics from the TRAINING set only
        mean_train <- colMeans(predictor_train$Z)
        sd_train <- apply(predictor_train$Z, 2, sd)
        
        # 2. Apply normalization to the TRAINING set
        predictor_train_normalized$Z <-
          scalar_normalize(
            predictor_train$Z,
            mean_train,
            sd_train
          )
        
        # 3. Apply normalization to the TEST set using TRAIN statistics
        predictor_test_normalized$Z <-
          scalar_normalize(
            predictor_test$Z,
            mean_train, # IMPORTANT: Use train mean
            sd_train    # IMPORTANT: Use train SD
          )
        
        return(
          list(
            predictor_train = predictor_train_normalized,
            predictor_test = predictor_test_normalized,
            mean_train = mean_train,
            sd_train = sd_train
          )
        )
      }
      ```
      
      
      ```{r}
      #' Scale the scalar predictors so that the variability between the
      #' functional and scalar predictors are comparable.
      #'
      #' This method calculates a scaling factor (sqrt(omega)) from the training
      #' set, defined as the ratio of the total sum of squares (L2 norms) of all
      #' functional components to the total sum of squares of all scalar components.
      #'
      #' @param train A predictor_hybrid object representing the training data.
      #' @param test A predictor_hybrid object representing the test data.
      #'
      #' @return A list containing the normalized 'train' and 'test' predictor_hybrid objects.
      #' @export
      btwn_normalize_train_test <- function(train, test){
        
        # Type checks
        stopifnot(
          "'train' must be of class 'predictor_hybrid'" = inherits(train, "predictor_hybrid"),
          "'test' must be of class 'predictor_hybrid'" = inherits(test, "predictor_hybrid")
        )
        
        # Create copies for modification
        train_normalized <- train
        test_normalized <- test
        
        # Z_s for clarity in calculation
        Z_s <- train$Z
        
        # --- 1. Calculate the numerator of omega (Sum of L2 norms of all functional samples in train) ---
        omega_numerator <- 0
        
        for (k in seq_len(train$n_functional)){ # for kth functional predictor
          # fda::inprod(fd, fd) returns the N x N matrix of inner products.
          # The diagonal elements are the L2 norm squared for each sample.
          functional_inprod_matrix <- fda::inprod(
            train$functional_list[[k]],
            train$functional_list[[k]]
          )
          
          # Sum the diagonal to get the total squared norm for this predictor k
          omega_numerator <- omega_numerator + sum(diag(functional_inprod_matrix))
        }
        
        # --- 2. Calculate the denominator of omega (Sum of squares of all scalar samples in train) ---
        omega_denominator <- sum(Z_s^2)
        
        # Handle case where scalar variance is zero (prevent division by zero)
        if (omega_denominator == 0) {
          warning("Scalar component sum of squares is zero. No scaling applied.")
          omega <- 1
        } else {
          # Calculate the scaling factor omega
          omega <- omega_numerator / omega_denominator
        }
        
        # --- 3. Apply scaling factor sqrt(omega) to both train and test scalar components ---
        scaling_factor <- sqrt(omega)
        
        train_normalized$Z <- scaling_factor * (train$Z)
        test_normalized$Z <- scaling_factor * (test$Z)
        
        return(
          list(
            predictor_train = train_normalized,
            predictor_test = test_normalized,
            scaling_factor = scaling_factor
          )
        )
      }
      
      ```
      
      
      ### normailize.all
      
      ```{r}
      #'
      #' @export
      split_and_normailize.all <- function(W_hybrid, response, train_ratio){
        
        
        # $W_hybrid =  W_hybrid_mock
        #$ response
        
        # train_ratio = TRAIN_RATIO
        
        # split.result <- split.all(W_hybrid, response, train_ratio)
        
        #     predictor_train = W_train,
        #   predictor_test = W_test ,
        #   response_train = response[train_idx],
        #   response_test = response[test_idx]
        
        # list(
        #   train = train_normalized,
        #   test = test_normalized,
        #   mean_train = mean_train,
        #   deno_train = deno_train
        #  ))
        
        curve_normalize_result <- curve_normalize_train_test(split.result$predictor_train, split.result$predictor_test)
        scalar_normalize_result <- scalar_normalize_train_test(curve_normalize_result$predictor_train, curve_normalize_result$predictor_test)
        btwn_normalize_result <- btwn_normalize_train_test(
          scalar_normalize_result$predictor_train,
          scalar_normalize_result$predictor_test
        )
        
        response_mean_train <-  mean(split.result$response_train)
        response_sd_train <-  sd(split.result$response_train)
        response_train <- (split.result$response_train - response_mean_train)/response_sd_train
        response_test <- (split.result$response_test - response_mean_train)/response_sd_train
        response_normalize_result = list(mean_train = response_mean_train, sd_train = response_sd_train)
        return(
          list(
            predictor_train = btwn_normalize_result$predictor_train,
            predictor_test = btwn_normalize_result$predictor_test,
            response_train = response_train,
            response_test = response_test,
            curve_normalize_result = curve_normalize_result,
            scalar_normalize_result = scalar_normalize_result,
            btwn_normalize_result = btwn_normalize_result,
            response_normalize_result = response_normalize_result
          )
        )
      }
      ```
      
      
      
      
      
      
      
      
      
      # Baseline method: FPCA regression
      
      
      ## Modified MFPCA Function for Simultaneous Training and Prediction
      
      The MFPCA (Multivariate Functional Principal Component Analysis) function computes multivariate functional PCA for functions defined across different (dimensional) domains.
      We have modified the original implementation to admit a new input, **`mFData_predict`**, allowing the function to calculate PCA scores for **both the training data and new test (or validation) data simultaneously** in a single run. This change was necessary because the standard `predict.MFPCAfit` function does not permit computing scores for new, out-of-sample functions.
      The multivariate functional principal component analysis inherently relies on a **univariate basis expansion** for each functional element ($X^{(j)}$). To simplify the code modification, we use only the univariate functional PCA option and have removed other options provided by the original MFPCA function. We modifies the main function and two internal functions.
      
      
      ### MFPCA
      
      ```{r}
      #'
      #' @export MFPCA
      MFPCA <- function(mFData, mFData_predict, M, uniExpansions, weights = rep(1, length(mFData)), fit = FALSE, approx.eigen = FALSE,
                        bootstrap = FALSE, nBootstrap = NULL, bootstrapAlpha = 0.05, bootstrapStrat = NULL,
                        verbose = options()$verbose)
      {
        if(! inherits(mFData, "multiFunData"))
          stop("Parameter 'mFData' must be passed as a multiFunData object.")
        
        # number of components
        p <- length(mFData)
        # number of observations
        N <- nObs(mFData)
        
        if(!all(is.numeric(M), length(M) == 1, M > 0))
          stop("Parameter 'M' must be passed as a number > 0.")
        
        if(!(is.list(uniExpansions) & length(uniExpansions) == p))
          stop("Parameter 'uniExpansions' must be passed as a list with the same length as 'mFData'.")
        
        if(!(is.numeric(weights) & length(weights) == p))
          stop("Parameter 'weights' must be passed as a vector with the same length as 'mFData'.")
        
        if(!is.logical(fit))
          stop("Parameter 'fit' must be passed as a logical.")
        
        if(!is.logical(approx.eigen))
          stop("Parameter 'approx.eigen' must be passed as a logical.")
        
        if(!is.logical(bootstrap))
          stop("Parameter 'bootstrap' must be passed as a logical.")
        
        if(bootstrap)
        {
          if(is.null(nBootstrap))
            stop("Specify number of bootstrap iterations.")
          
          if(any(!(0 < bootstrapAlpha & bootstrapAlpha < 1)))
            stop("Significance level for bootstrap confidence bands must be in (0,1).")
          
          if(!is.null(bootstrapStrat))
          {
            if(!is.factor(bootstrapStrat))
              stop("bootstrapStrat must be either NULL or a factor.")
            
            if(length(bootstrapStrat) != nObs(mFData))
              stop("bootstrapStrat must have the same length as the number of observations in the mFData object.")
          }
        }
        
        if(!is.logical(verbose))
          stop("Parameter 'verbose' must be passed as a logical.")
        
        # dimension for each component
        dimSupp <- dimSupp(mFData)
        
        # get type of univariate expansions
        type <- vapply(uniExpansions, function(l){l$type}, FUN.VALUE = "")
        
        # de-mean functions -> coefficients are also de-meaned!
        # do not de-mean in uFPCA, as PACE gives a smooth estimate of the mean (see below)
        m <- meanFunction(mFData, na.rm = TRUE) # ignore NAs in data
        for(j in seq_len(p))
        {
          if(type[j] != "uFPCA")
            mFData[[j]] <- mFData[[j]] - m[[j]]
        }
        
        if(verbose)
          cat("Calculating univariate basis expansions (", format(Sys.time(), "%T"), ")\n", sep = "")
        
        
        
        ########### Below is modified by Jongmin Mun #############################
        # calculate univariate basis expansion for all components
        uniBasis <- mapply(
          function(expansion, data, data_predict){
            do.call(univDecomp, c(
              list(funDataObject = data, funDataObject_predict = data_predict),
              expansion)
            )
          },
          expansion = uniExpansions,
          data = mFData,
          data_predict = mFData_predict,
          SIMPLIFY = FALSE
        )
        ########### Above is modified by Jongmin Mun #############################
        
        
        
        
        
        # for uFPCA: replace estimated mean in m
        for(j in seq_len(p))
        {
          if(type[j] == "uFPCA")
            m[[j]] <- uniBasis[[j]]$meanFunction
        }
        
        # Multivariate FPCA
        npc <- vapply(uniBasis, function(x){dim(x$scores)[2]}, FUN.VALUE = 0) # get number of univariate basis functions
        
        if(M > sum(npc))
        {
          M <- sum(npc)
          warning("Function MFPCA: total number of univariate basis functions is smaller than given M. M was set to ", sum(npc), ".")
        }
        
        # check if non-orthonormal basis functions used
        if(all(foreach::foreach(j = seq_len(p), .combine = "c")%do%{uniBasis[[j]]$ortho}))
          Bchol = NULL
        else
        {
          # Cholesky decomposition of B = block diagonal of Cholesky decompositions
          Bchol <- Matrix::bdiag(lapply(uniBasis, function(l){
            if(l$ortho)
              res <- Matrix::Diagonal(n = ncol(l$scores))
            else
              res <- Matrix::chol(l$B)
            
            return(res)}))
        }
        
        if(verbose)
          cat("Calculating MFPCA (", format(Sys.time(), "%T"), ")\n", sep = "")
        
        mArgvals <- if (utils::packageVersion("funData") <= "1.2") {
          getArgvals(mFData)
        } else {
          funData::argvals(mFData)
        }
        
        res <- calcMFPCA(N = N, p = p, Bchol = Bchol, M = M, type = type, weights = weights,
                         npc = npc, argvals = mArgvals, uniBasis = uniBasis, fit = fit, approx.eigen = approx.eigen)
        
        res$meanFunction <- m # return mean function, too
        
        names(res$functions) <- names(mFData)
        
        #############modified by jongmin mun
        
        
        #####################################
        if(fit)
        {
          res$fit <- m + res$fit # add mean function to fits
          names(res$fit) <- names(mFData)
        }
        
        # give correct names
        namesList <- lapply(mFData, names)
        if(!all(vapply(namesList, FUN = is.null, FUN.VALUE = TRUE)))
        {
          if(length(unique(namesList)) != 1)
            warning("Elements have different curve names. Use names of the first element for the results.")
          
          row.names(res$scores) <- namesList[[1]]
          
          if(fit)
            for(i in seq_len(p))
              names(res$fit[[i]]) <- namesList[[1]]
        }
        
        
        if(type[1] == "uFPCA"){
          scores.pred<-calcMFPCA_predict(N = nObs(mFData_predict),
                                         p = p,
                                         M = M,
                                         weights = weights,
                                         npc = npc,
                                         uniBasis = uniBasis,
                                         normFactors = res$normFactors,
                                         vectors_train=res$vectors,
                                         values_train=res$values
          )
        }
        res$scores.pred <- scores.pred
        
        class(res) <- "MFPCAfit"
        return(res)
      }
      
      ```
      
      
      ### fpcaBasis
      
      ```{r}
      ################# modified by Jongmin Mun #####################
      fpcaBasis <- function(funDataObject, funDataObject_predict, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
      {
        FPCA <- PACE(funDataObject, predData = NULL, nbasis, pve, npc, makePD, cov.weight.type)
        FPCA_pred <- PACE(funDataObject = funDataObject, predData = funDataObject_predict, nbasis, pve, npc, makePD, cov.weight.type)
        return(list(scores = FPCA$scores,
                    scores.pred = FPCA_pred$scores,
                    ortho = TRUE,
                    functions = FPCA$functions,
                    meanFunction = FPCA$mu
        ))
      }
      ```
      
      
      ### univDecomp
      ```{r}
      univDecomp <- function(type, funDataObject, funDataObject_predict, ...)
      {
        # Parameter checking
        if(is.null(type))
          stop("Parameter 'type' is missing.")
        else
        {
          if(!is.character(type))
            stop("Parameter 'type' must be a character string. See ?univDecomp for details.")
        }
        
        if(class(funDataObject) != "funData")
          stop("Parameter 'funDataObject' must be a funData object.")
        
        
        # get all arguments (except for function call and type)
        params <- list(...)
        
        # check if type and data are of correct type
        if(is.null(type))
          stop("univDecomp: must specify 'type'.")
        
        if(!inherits(type, "character"))
          stop("univDecomp: 'type' must be of class character.")
        
        if(is.null(funDataObject))
          stop("univDecomp: must specify 'funDataObject'.")
        
        if(class(funDataObject) != "funData")
          stop("univDecomp: 'funDataObject' must be of class funData.")
        
        params$funDataObject <- funDataObject # add funDataObject (-> make sure is evaluated in correct env.)
        
        ############### modified by Jongmin Mun ######################
        if (type == "uFPCA"){params$funDataObject_predict <- funDataObject_predict}
        ############### modified by Jongmin Mun ######################
        
        res <- switch(type,
                      "given" = do.call(givenBasis, params),
                      "uFPCA" = do.call(fpcaBasis, params),
                      "UMPCA" = do.call(umpcaBasis, params),
                      "FCP_TPA" = do.call(fcptpaBasis, params),
                      "splines1D" = do.call(splineBasis1D, params),
                      "splines1Dpen" = do.call(splineBasis1Dpen, params),
                      "splines2D" = do.call(splineBasis2D, params),
                      "splines2Dpen" = do.call(splineBasis2Dpen, params),
                      "fda" = do.call(fdaBasis, params),
                      "DCT2D" = do.call(dctBasis2D, params),
                      "DCT3D" = do.call(dctBasis3D, params),
                      stop("Univariate Decomposition for 'type' = ", type, " not defined!")
        )
        
        if(res$ortho == FALSE & is.null(res$B))
          stop("UnivDecomp: must provide integral matrix B for non-orthonormal basis functions.")
        
        return(res)
      }
      ```
      
      
      
      ## Unmodified functions
      
      ### calcBasisIntegrals
      ```{r}
      
      # define global variable j, used by the foreach package and confusing R CMD CHECK
      globalVariables('j')
      
      #' Utility function that calculates matrix of basis-scalar products (one dimension)
      #'
      #' If the element \eqn{X^{(j)}}{X^{(j)}} is expanded in basis functions \eqn{b_i^{(j)}(t),~ i = 1, \ldots, K_j}{b_i(t)},
      #' this function calculates the \eqn{K_j \times K_j}{K_j  \times K_j} matrix \eqn{B^{(jj)}}{B^{(jj)}} with entries
      #' \deqn{B^{(jj)}_{mn} = \int_{\mathcal{T_j}} b_m^{(j)}(t) b_n^{(j)}(t) \mathrm{d} t}.
      #'
      #' @section Warning: This function is implemented only for functions on one- or two-dimensional domains.
      #'
      #' @param basisFunctions Array of \code{npc} basis functions of dimensions \code{npc x M1} or \code{npc x M1 x M2}.
      #' @param dimSupp dimension of the support of the basis functions (1 or 2)
      #' @param argvals List of corresponding x-values.
      #'
      #' @return A matrix containing the scalar product of all combinations of basis functions (matrix \eqn{B^{(j)}})
      #'
      #' @seealso \code{\link{MFPCA}}, \code{\link[funData]{dimSupp}}
      #'
      #' @keywords internal
      calcBasisIntegrals <- function(basisFunctions, dimSupp, argvals)
      {
        npc <- dim(basisFunctions)[1]
        
        #  integral basis matrix
        B <- array(0, dim = c(npc, npc))
        
        
        if(dimSupp == 1) # one-dimensional domain
        {
          w <- funData::.intWeights(argvals[[1]])
          
          for(m in seq_len(npc))
          {
            for(n in seq_len(m))
              B[m, n] <- B[n, m] <- (basisFunctions[m, ]* basisFunctions[n, ])%*%w
          }
        }
        else # two-dimesional domain (otherwise function stops before!)
        {
          w1 <- t(funData::.intWeights(argvals[[1]]))
          w2 <- funData::.intWeights(argvals[[2]])
          
          for(m in seq_len(npc))
          {
            for(n in seq_len(m))
              B[m, n] <- B[n, m] <-  w1 %*%(basisFunctions[m, , ]* basisFunctions[n, ,])%*%w2
          }
        }
        
        return(B)
      }
      ```
      
      
      
      
      ### calcMFPCA
      ```{r}
      #' Internal function that implements the MFPCA algorithm for given univariate decompositions
      #' @export
      calcMFPCA <- function(N, p, Bchol, M, type, weights, npc, argvals, uniBasis, fit = FALSE, approx.eigen = FALSE)
      {
        # combine all scores
        allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores}
        
        # block vector of weights
        allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}
        
        Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)
        
        # check if approximation is appropriate (cf. irlba)
        if(approx.eigen & (M > min(N, sum(npc))/2))
        {
          warning("Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
          approx.eigen = FALSE
        }
        
        # check if non-orthonormal basis functions used and calculate PCA on scores
        if(is.null(Bchol))
        {
          if(approx.eigen)
          {
            tmpSVD <- irlba::irlba(as.matrix(Z), nv = M)
            
            vectors <- tmpSVD$v
            values <- tmpSVD$d[seq_len(M)]^2
          }
          else
          {
            if(sum(npc) > 1000)
              warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")
            
            e <- eigen(stats::cov(allScores) * outer(allWeights, allWeights, "*"))
            
            values <- e$values[seq_len(M)]
            vectors <- e$vectors[,seq_len(M)]
          }
        }
        else
        {
          if(approx.eigen)
          {
            tmpSVD <- irlba::irlba(as.matrix(Matrix::tcrossprod(Z, Bchol)), nv = M)
            
            vectors <- Matrix::crossprod(Bchol, tmpSVD$v)
            values <- tmpSVD$d[seq_len(M)]^2
          }
          else
          {
            if(sum(npc) > 1000)
              warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")
            
            e <- eigen(Matrix::crossprod(Bchol) %*% (stats::cov(allScores) * outer(allWeights, allWeights, "*")))
            
            values <- Re(e$values[seq_len(M)])
            vectors <- Re(e$vectors[,seq_len(M)])
          }
        }
        
        # normalization factors
        normFactors <- 1/sqrt(diag(as.matrix(Matrix::crossprod(Z %*% vectors))))
        
        ### Calculate scores
        scores <- Z %*% vectors * sqrt(N-1) # see defintion of Z above!
        scores <- as.matrix(scores %*% diag(sqrt(values) * normFactors, nrow = M, ncol = M)) # normalization
        
        ### Calculate eigenfunctions (incl. normalization)
        npcCum <- cumsum(c(0, npc)) # indices for blocks (-1)
        
        tmpWeights <- as.matrix(Matrix::crossprod(Z, Z %*%vectors))
        eFunctions <- foreach::foreach(j = seq_len(p)) %do% {
          univExpansion(type = type[j],
                        scores = 1/sqrt(weights[j] * values) * normFactors * t(tmpWeights[npcCum[j]+seq_len(npc[j]), , drop = FALSE]),
                        argvals = argvals[[j]],
                        functions = uniBasis[[j]]$functions,
                        params = uniBasis[[j]]$settings)
        }
        
        res <- list(values = values,
                    functions = multiFunData(eFunctions),
                    scores = scores,
                    vectors = vectors,
                    values = values,
                    normFactors = normFactors,
                    uniBasis = uniBasis,
                    uniExpansions = uniExpansions
        )
        
        # calculate truncated Karhunen-Loeve representation (no mean here)
        if(fit)
          res$fit <- multivExpansion(multiFuns = res$functions, scores = scores)
        
        return(res)
      }
      ```
      
      
      ### calcMFPCA
      
      ```{r}
      #' Internal function that implements the MFPCA algorithm for given univariate decompositions
      #' @export
      #'
      calcMFPCA <- function(N, p, Bchol, M, type, weights, npc, argvals, uniBasis, fit = FALSE, approx.eigen = FALSE)
      {
        # combine all scores
        allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores}
        
        # block vector of weights
        allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}
        
        Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)
        
        # check if approximation is appropriate (cf. irlba)
        if(approx.eigen & (M > min(N, sum(npc))/2))
        {
          warning("Calculating a large percentage of principal components, approximation may not be appropriate.
            'approx.eigen' set to FALSE.")
          approx.eigen = FALSE
        }
        
        # check if non-orthonormal basis functions used and calculate PCA on scores
        if(is.null(Bchol))
        {
          if(approx.eigen)
          {
            tmpSVD <- irlba::irlba(as.matrix(Z), nv = M)
            
            vectors <- tmpSVD$v
            values <- tmpSVD$d[seq_len(M)]^2
          }
          else
          {
            if(sum(npc) > 1000)
              warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")
            
            e <- eigen(stats::cov(allScores) * outer(allWeights, allWeights, "*"))
            
            values <- e$values[seq_len(M)]
            vectors <- e$vectors[,seq_len(M)]
          }
        }
        else
        {
          if(approx.eigen)
          {
            tmpSVD <- irlba::irlba(as.matrix(Matrix::tcrossprod(Z, Bchol)), nv = M)
            
            vectors <- Matrix::crossprod(Bchol, tmpSVD$v)
            values <- tmpSVD$d[seq_len(M)]^2
          }
          else
          {
            if(sum(npc) > 1000)
              warning("MFPCA with > 1000 univariate eigenfunctions and approx.eigen = FALSE. This may take some time...")
            
            e <- eigen(Matrix::crossprod(Bchol) %*% (stats::cov(allScores) * outer(allWeights, allWeights, "*")))
            
            values <- Re(e$values[seq_len(M)])
            vectors <- Re(e$vectors[,seq_len(M)])
          }
        }
        
        # normalization factors
        normFactors <- 1/sqrt(diag(as.matrix(Matrix::crossprod(Z %*% vectors))))
        
        ### Calculate scores
        scores <- Z %*% vectors * sqrt(N-1) # see defintion of Z above!
        scores <- as.matrix(scores %*% diag(sqrt(values) * normFactors, nrow = M, ncol = M)) # normalization
        
        ### Calculate eigenfunctions (incl. normalization)
        npcCum <- cumsum(c(0, npc)) # indices for blocks (-1)
        
        tmpWeights <- as.matrix(Matrix::crossprod(Z, Z %*%vectors))
        eFunctions <- foreach::foreach(j = seq_len(p)) %do% {
          univExpansion(type = type[j],
                        scores = 1/sqrt(weights[j] * values) * normFactors * t(tmpWeights[npcCum[j]+seq_len(npc[j]), , drop = FALSE]),
                        argvals = argvals[[j]],
                        functions = uniBasis[[j]]$functions,
                        params = uniBasis[[j]]$settings)
        }
        
        res <- list(values = values,
                    functions = multiFunData(eFunctions),
                    scores = scores,
                    vectors = vectors,
                    values = values,
                    normFactors = normFactors,
                    uniBasis = uniBasis,
                    uniExpansions = uniExpansions
        )
        
        # calculate truncated Karhunen-Loeve representation (no mean here)
        if(fit)
          res$fit <- multivExpansion(multiFuns = res$functions, scores = scores)
        
        return(res)
      }
      ```
      
      ### calcMFPCA_predict
      
      ```{r}
      #' calcMFPCA_predict
      #' @export
      #'
      calcMFPCA_predict <- function(N, p, M, weights, npc, uniBasis,
                                    normFactors_train, vectors_train, values_train)
      {
        # combine all scores
        allScores <- foreach::foreach(j = seq_len(p), .combine = "cbind")%do%{uniBasis[[j]]$scores.pred}
        
        # block vector of weights
        allWeights <- foreach::foreach(j = seq_len(p), .combine = "c")%do%{rep(sqrt(weights[j]), npc[j])}
        
        Z <- allScores %*% Matrix::Diagonal(x = allWeights) / sqrt(N-1)
        
        
        ### Calculate scores
        scores <- Z %*% vectors_train * sqrt(N-1) # see defintion of Z above!
        scores <- as.matrix(scores %*% diag(sqrt(values_train) * normFactors_train, nrow = M, ncol = M)) # normalization
        
        return(scores)
      }
      
      ```
      
      
      ### .PACE
      
      ```{r}
      .PACE <- function(X, Y, Y.pred = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
      {
        if (is.null(Y.pred))
          Y.pred = Y
        D = NCOL(Y)
        if(D != length(X)) # check if number of observation points in X & Y are identical
          stop("different number of (potential) observation points differs in X and Y!")
        I = NROW(Y)
        I.pred = NROW(Y.pred)
        d.vec = rep(X, each = I) # use given X-values for estimation of mu
        gam0 = mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
        mu = mgcv::predict.gam(gam0, newdata = data.frame(d.vec = X))
        Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)
        cov.sum = cov.count = cov.mean = matrix(0, D, D)
        for (i in seq_len(I)) {
          obs.points = which(!is.na(Y[i, ]))
          cov.count[obs.points, obs.points] = cov.count[obs.points, obs.points] + 1
          cov.sum[obs.points, obs.points] = cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
        }
        G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)
        diag.G0 = diag(G.0)
        diag(G.0) = NA
        row.vec = rep(X, each = D) # use given X-values
        col.vec = rep(X, D) # use given X-values
        cov.weights <- switch(cov.weight.type,
                              none = rep(1, D^2),
                              counts = as.vector(cov.count),
                              stop("cov.weight.type ", cov.weight.type, " unknown in smooth covariance estimation"))
        
        npc.0 = matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0)~te(row.vec, col.vec, k = nbasis), weights = cov.weights),
                                         newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
        npc.0 = (npc.0 + t(npc.0))/2
        # no extra-option (useSymm) as in fpca.sc-method
        if (makePD) { # see fpca.sc
          npc.0 <- {
            tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE,
                                  do2eigen = TRUE, trace = options()$verbose)
            as.matrix(tmp$mat)
          }
        }
        
        ### numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Chapter 8)
        w <- funData::.intWeights(X, method = "trapezoidal")
        Wsqrt <- diag(sqrt(w))
        Winvsqrt <- diag(1/(sqrt(w)))
        V <- Wsqrt %*% npc.0 %*% Wsqrt
        evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values
        ###
        evalues = replace(evalues, which(evalues <= 0), 0)
        npc = ifelse(is.null(npc), min(which(cumsum(evalues)/sum(evalues) > pve)), npc)
        efunctions = matrix(Winvsqrt%*%eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = D, ncol = npc)
        evalues = eigen(V, symmetric = TRUE, only.values = TRUE)$values[seq_len(npc)]  # use correct matrix for eigenvalue problem
        cov.hat = efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
        ### numerical integration for estimation of sigma2
        T.len <- X[D] - X[1] # total interval length
        T1.min <- min(which(X >= X[1] + 0.25*T.len)) # left bound of narrower interval T1
        T1.max <- max(which(X <= X[D] - 0.25*T.len)) # right bound of narrower interval T1
        DIAG = (diag.G0 - diag(cov.hat))[T1.min :T1.max] # function values
        # weights
        w <- funData::.intWeights(X[T1.min:T1.max], method = "trapezoidal")
        sigma2 <- max(1/(X[T1.max]-X[T1.min]) * sum(DIAG*w, na.rm = TRUE), 0) #max(1/T.len * sum(DIAG*w), 0)
        ####
        D.inv = diag(1/evalues, nrow = npc, ncol = npc)
        Z = efunctions
        Y.tilde = Y.pred - matrix(mu, I.pred, D, byrow = TRUE)
        fit = matrix(0, nrow = I.pred, ncol = D)
        scores = matrix(NA, nrow = I.pred, ncol = npc)
        # no calculation of confidence bands, no variance matrix
        for (i.subj in seq_len(I.pred)) {
          obs.points = which(!is.na(Y.pred[i.subj, ]))
          if (sigma2 == 0 & length(obs.points) < npc) {
            stop("Measurement error estimated to be zero and there are fewer observed points than PCs; scores cannot be estimated.")
          }
          Zcur = matrix(Z[obs.points, ], nrow = length(obs.points),
                        ncol = dim(Z)[2])
          ZtZ_sD.inv = solve(crossprod(Zcur) + sigma2 * D.inv)
          scores[i.subj, ] = ZtZ_sD.inv %*% crossprod(Zcur, Y.tilde[i.subj, obs.points])
          fit[i.subj, ] = t(as.matrix(mu)) + tcrossprod(scores[i.subj, ], efunctions)
        }
        ret.objects = c("fit", "scores", "mu", "efunctions", "evalues",
                        "npc", "sigma2") # add sigma2 to output
        ret = lapply(seq_len(length(ret.objects)), function(u) get(ret.objects[u]))
        names(ret) = ret.objects
        ret$estVar <- diag(cov.hat)
        return(ret)
      }
      ```
      
      ### PACE
      
      ```{r}
      PACE <- function(funDataObject, predData = NULL, nbasis = 10, pve = 0.99, npc = NULL, makePD = FALSE, cov.weight.type = "none")
      {
        # check inputs
        if(! class(funDataObject) %in% c("funData", "irregFunData"))
          stop("Parameter 'funDataObject' must be a funData or irregFunData object.")
        if(dimSupp(funDataObject) != 1)
          stop("PACE: Implemented only for funData objects with one-dimensional support.")
        if(methods::is(funDataObject, "irregFunData")) # for irregular functional data, use funData representation
          funDataObject <- as.funData(funDataObject)
        
        if(is.null(predData))
          Y.pred = NULL # use only funDataObject
        else
        {
          if(!isTRUE(all.equal(funDataObject@argvals, predData@argvals)))
            stop("PACE: funDataObject and predData must be defined on the same domains!")
          
          Y.pred = predData@X
        }
        
        
        if(!all(is.numeric(nbasis), length(nbasis) == 1, nbasis > 0))
          stop("Parameter 'nbasis' must be passed as a number > 0.")
        
        if(!all(is.numeric(pve), length(pve) == 1, 0 <= pve, pve <= 1))
          stop("Parameter 'pve' must be passed as a number between 0 and 1.")
        
        if(!is.null(npc) & !all(is.numeric(npc), length(npc) == 1, npc > 0))
          stop("Parameter 'npc' must be either NULL or passed as a number > 0.")
        
        if(!is.logical(makePD))
          stop("Parameter 'makePD' must be passed as a logical.")
        
        if(!is.character(cov.weight.type))
          stop("Parameter 'cov.weight.type' must be passed as a character.")
        
        
        res <- .PACE(X = funDataObject@argvals[[1]],
                     Y = funDataObject@X,
                     Y.pred = Y.pred,
                     nbasis = nbasis, pve = pve, npc = npc, makePD = makePD,
                     cov.weight.type = cov.weight.type)
        
        return(list(mu = funData(funDataObject@argvals, matrix(res$mu, nrow = 1)),
                    values = res$evalues,
                    functions = funData(funDataObject@argvals, t(res$efunctions)),
                    scores = res$scores,
                    scores.pred = res$scores.pred,
                    fit = funData(funDataObject@argvals, res$fit),
                    npc = res$npc,
                    sigma2 = res$sigma2,
                    estVar = funData(funDataObject@argvals, matrix(res$estVar, nrow = 1))
        ))
      }
      ```
      
      ### expandBasisFunction
      
      ```{r}
      expandBasisFunction <- function(scores, argvals = functions@argvals, functions)
      {
        if(dim(scores)[2] != nObs(functions))
          stop("expandBasisFunction: number of scores for each observation and number of eigenfunctions does not match.")
        
        # collapse higher-dimensional functions, multiply with scores and resize the result
        d <- dim(functions@X)
        nd <- length(d)
        
        if(nd == 2)
          resX <- scores %*% functions@X
        
        if(nd == 3)
        {
          resX <- array(NA, dim = c(dim(scores)[1], d[-1]))
          
          for(i in seq_len(d[2]))
            resX[,i,] <- scores %*% functions@X[,i,]
        }
        
        if(nd == 4)
        {
          resX <- array(NA, dim = c(dim(scores)[1], d[-1]))
          
          for(i in seq_len(d[2]))
            for(j in seq_len(d[3]))
              resX[,i,j,] <- scores %*% functions@X[,i,j,]
        }
        
        if(nd > 4) # slow solution due to aperm
        {
          resX <- aperm(plyr::aaply(.data = functions@X, .margins = 3:nd,
                                    .fun = function(x,y){y %*% x}, y = scores),
                        c(nd-1,nd, seq_len((nd-2))))
          dimnames(resX) <- NULL
        }
        
        return( funData(argvals, resX) )
      }
      ```
      
      
      ### univExpansion
      
      ```{r}
      univExpansion <- function(type, scores, argvals = ifelse(!is.null(functions), functions@argvals, NULL), functions, params = NULL)
      {
        # Parameter checking
        if(is.null(type))
          stop("Parameter 'type' is missing.")
        else
        {
          if(!is.character(type))
            stop("Parameter 'type' must be a character string. See ?univExpansion for details.")
        }
        
        if(is.null(scores))
          stop("Parameter 'scores' is missing.")
        else
        {
          if(!is.matrix(scores))
            stop("Parameter 'scores' must be passed as a matrix.")
        }
        
        if(is.numeric(argvals))
        {
          argvals <- list(argvals)
          warning("Parameter 'argvals' was passed as a vector and transformed to a list.")
        }
        
        if(is.null(functions))
        {
          if(is.null(argvals))
            stop("Must pass 'argvals' if 'functions' is NULL.")
          else
          {
            if(!is.list(argvals))
              stop("Parameter 'argvals' must be passed as a list.")
          }
        }
        else
        {
          if(class(functions) != "funData")
            stop("Parameter 'functions' must be a funData object.")
          
          # check interaction with other parameters
          if(nObs(functions) != NCOL(scores))
            stop("Number of scores per curve does not match the number of basis functions.")
          
          if(!is.null(argvals) & !isTRUE(all.equal(argvals, functions@argvals)))
            stop("The parameter 'argvals' does not match the argument values of 'functions'.")
        }
        
        if(!is.null(params) & !is.list(params))
          stop("The parameter 'params' must be passed as a list.")
        
        # start calculations
        params$scores <- scores
        params$functions <- functions
        
        if(is.numeric(argvals))
          argvals <- list(argvals)
        
        params$argvals <- argvals
        
        res <- switch(type,
                      "given" = do.call(expandBasisFunction, params),
                      "uFPCA" = do.call(expandBasisFunction, params),
                      "UMPCA" = do.call(expandBasisFunction, params),
                      "FCP_TPA" = do.call(expandBasisFunction, params),
                      "splines1D" = do.call(splineFunction1D, params),
                      "splines1Dpen" = do.call(splineFunction1D, params),
                      "splines2D" = do.call(splineFunction2D, params),
                      "splines2Dpen" = do.call(splineFunction2Dpen, params),
                      "fda" = do.call(expandBasisFunction, params),
                      "DCT2D" = do.call(dctFunction2D, params),
                      "DCT3D" = do.call(dctFunction3D, params),
                      "default" = do.call(expandBasisFunction, params),
                      stop("Univariate Expansion for 'type' = ", type, " not defined!")
        )
        
        return(res)
      }
      ```
      
      ## fit.mfpca_regression
      Performs sample splitting, normalization, and fits models for every combination of functional and scalar principal components up to n_pc_max. Selects the best model based on test set RMSE.
      
      ```{r}
      # --- 1. MAIN FUNCTION: fit.mfpca_regression ---
      
      #' End-to-End MFPCA Regression Pipeline
      #'
      #' Performs sample splitting, normalization, cross-validation (CV) to select
      #' optimal functional (FPC) and scalar (SPC) components, linear regression fitting,
      #' and test set prediction.
      #'
      #' @param W predictor_hybrid object containing all predictors.
      #' @param y Vector. Response variable.
      #' @param train_ratio Numeric scalar between 0 and 1 for splitting data (default 0.7).
      #' @param n_pc_max Integer. Max number of FPCs/PCs to consider (must be <= min(nbasis, n_scalar)).
      #' @return List containing the best model fit, optimal component count, and a list of all test RMSEs.
      #' @export
      fit.mfpca_regression <- function(W, y, train_ratio = 0.7, n_pc_max = 10) {
        
        # --- 0. Setup and Initial Checks ---
        stopifnot(inherits(W, "predictor_hybrid"))
        
        # Extract setup parameters
        eval_point <- W$eval_point
        n_basis <- W$n_basis_list[[1]]
        n_scalar <- W$n_scalar
        
        # Adjust max components to safe bounds
        max_comp <- min(n_pc_max, n_basis, n_scalar)
        if (n_pc_max != max_comp) {
          warning(paste("n_pc_max adjusted from", n_pc_max, "to", max_comp,
                        "due to limits in nbasis or n_scalar."))
          n_pc_max <- max_comp
        }
        
        if (n_pc_max < 1) stop("n_pc_max must be at least 1.")
        
        uniExpansions <- list(
          list(type = "uFPCA", nbasis = n_basis, npc = n_pc_max),
          list(type = "uFPCA", nbasis = n_basis, npc = n_pc_max)
        )
        
        
        # --- 1. Split and Normalize Full Data ---
        split_normalize_result <- split_and_normailize.all(W, y, train_ratio)
        f1_train <- split_normalize_result$predictor_train$functional_list[[1]]
        f2_train <- split_normalize_result$predictor_train$functional_list[[2]]
        f1_test <- split_normalize_result$predictor_test$functional_list[[1]]
        f2_test <- split_normalize_result$predictor_test$functional_list[[2]]
        
        
        # Format functional data
        predictor_train <- multiFunData(fd2funData(f1_train, eval_point), fd2funData(f2_train, eval_point))
        predictor_test <- multiFunData(fd2funData(f1_test, eval_point), fd2funData(f2_test, eval_point))
        
        y_train_norm <- split_normalize_result$response_train
        Z_train_norm <- split_normalize_result$predictor_train$Z
        Z_test_norm <- split_normalize_result$predictor_test$Z
        y_test_norm <- split_normalize_result$response_test
        
        # Fit MFPCA (M=n_pc_max to get max scores)
        MFPCA_fit <- MFPCA(mFData = predictor_train, mFData_predict = predictor_test,
                           M = n_pc_max, uniExpansions = uniExpansions)
        
        # Fit Scalar PCA (prcomp)
        PCA_fit <- prcomp(Z_train_norm, center = FALSE, scale = FALSE)
        
        # Extract max scores
        FPC_score_train_max <- MFPCA_fit$scores
        FPC_score_test_max <- MFPCA_fit$scores.pred
        PC_score_train_max <- PCA_fit$x
        PC_score_test_max <- predict(PCA_fit, Z_test_norm)
        
        
        # --- 3. Iterate and Fit Models (Testing FPC=npc, PC=npc) ---
        
        linear_fit_list <- list()
        test_rmse_list <- numeric(n_pc_max)
        
        for (n_pc in 1:n_pc_max) {
          
          # 3a. Select required number of components
          FPC_train_select <- FPC_score_train_max[, 1:n_pc, drop = FALSE]
          FPC_test_select <- FPC_score_test_max[, 1:n_pc, drop = FALSE]
          PC_train_select <- PC_score_train_max[, 1:n_pc, drop = FALSE]
          PC_test_select <- PC_score_test_max[, 1:n_pc, drop = FALSE]
          
          # Create column names
          col_names_fpc <- paste0("FPC", 1:n_pc)
          col_names_pc <- paste0("PC", 1:n_pc)
          
          # 3b. Prepare data frames for linear model
          train_data_lm <- data.frame(
            response = y_train_norm,
            FPC_train_select,
            PC_train_select
          )
          colnames(train_data_lm)[-1] <- c(col_names_fpc, col_names_pc)
          
          test_data_lm <- data.frame(
            FPC_test_select,
            PC_test_select
          )
          colnames(test_data_lm) <- c(col_names_fpc, col_names_pc)
          
          # 3c. Fit Linear Regression
          linear_fit_list[[n_pc]] <- lm(response ~ ., data = train_data_lm)
          
          # 3d. Predict and Calculate RMSE on Test Set
          y_pred_test <- predict(linear_fit_list[[n_pc]], test_data_lm)
          test_rmse_list[n_pc] <- sqrt(mean((y_pred_test - y_test_norm)^2))
        }
        
        # --- 4. Select Best Model ---
        n_pc_best <- which.min(test_rmse_list)
        best_rmse <- test_rmse_list[n_pc_best]
        best_fit <- linear_fit_list[[n_pc_best]]
        
        # --- 5. Return Results ---
        return(list(
          n_pc_best = n_pc_best,
          best_test_rmse = best_rmse,
          test_rmse_by_n_pc = test_rmse_list,
          final_model = best_fit,
          model_list = linear_fit_list,
          regression_coefficients = coef(best_fit)
        ))
      }
      ```
      
      
      # data generation
      
      ## matrix normal
      ```{r}
      trial_name <- "matrix_normal"
      
      setwd("~")
      if (
        substr(getwd(), nchar(getwd())-10+1, nchar(getwd())) == "/Documents"
      ){
        wd_now <-  substr(getwd(), 1, nchar(getwd())-10)
      }else{
        wd_now <- getwd()
      }
      packagepath <- paste0( wd_now, "/Documents/GitHub/fsPLS") #mac
      datapath <- paste(
        wd_now,
        paste("Documents/GitHub/fsPLS/simul/data_generation", trial_name, sep = "/"),
        sep = "/"
      )
      
      devtools::load_all(packagepath)
      
      ##### parameter settings by user ######
      n_rep    <- 100
      n_sample <- 100
      n_eval            <- 99
      n_predictor_scalar <- 16
      test.ratio <- 0.3
      train.ratio <- 1 - test.ratio
      # parameters for predictor generation
      ar_slope <- 0.9
      gaussian_bandwidth <- 15
      freq_1 <- 10
      freq_2 <- 20
      signal_noise_ratio <- 0.9
      scalar_signal_ratio <- 0.1
      
      n_basis_predictor <- 20
      
      # parameters for coefficient generation
      n_basis_coef      <- 30
      cycle_1 <- 5
      cycle_2 <- 1
      
      reg_coef_1_formula <- function(x){2 * sin(0.5*cycle_1*pi*x) + 4* sin(1.5*cycle_1*pi*x) + 5 * sin( 2.5 * cycle_1*pi * x)}
      reg_coef_2_formula <- function(x){
        exp(
          -(cycle_2*(x-0.2)^2)/(0.04^2) )-
          exp( -(cycle_2*(x-0.4)^2)/(0.04^2) ) +
          exp( -(cycle_2*(x-0.6)^2)/(0.04^2) ) -
          exp(  -(cycle_2*(x-0.8)^2)/(0.04^2) )}
      
      
      
      ##########################################
      ##########################################
      
      # predictor generating parameters
      eval_point      <- seq(0,1, length = n_eval)
      basis_predictor <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_predictor)
      lambda_pls_base <- max(abs(inprod(basis_predictor, basis_predictor)))/max(abs(getbasispenalty(basis_predictor)))
      
      ## functional regression coefficient
      basis_reg_coef  <- create.bspline.basis(rangeval = c(0,1), nbasis = n_basis_coef)
      reg_coef_1 <- fda::Data2fd(argvals = eval_point, y = as.vector(reg_coef_1_formula(eval_point)), basisobj = basis_reg_coef)
      reg_coef_2 <- fda::Data2fd(argvals = eval_point, y = as.vector(reg_coef_2_formula(eval_point)), basisobj = basis_reg_coef)
      #
      ## scalar regression coefficient
      set.seed(1); scalar_reg_coef <- round( rnorm(n_predictor_scalar, 0, 1), 1 )
      
      
      
      
      
      
      
      # part II. data generation
      data_fun_1 <- data_fun_2<- data_scalar <- signal_functional <- signal_scalar <- signal_total <- response <- list()
      U_half <- matrix(rnorm(9), nrow=3)
      #U = U_half %*% t(U_half)
      #U
      U <- matrix(
        c(0.3585408, 0.3821649, 0.4974333,
          0.3821649, 5.6928664, 0.2608125,
          0.4974333, 0.2608125, 1.6703254
        ), nrow = 3)
      for (rep_number in 1:n_rep){
        set.seed(rep_number*10)
        vec_1 <- matrix(NA, nrow=n_sample, ncol = n_eval)
        vec_2 <- matrix(NA, nrow=n_sample, ncol = n_eval)
        vec_3 <- matrix(NA, nrow=n_sample, ncol = n_eval)
        for (obs in 1:n_sample){
          realization <- rmatrixnorm(1,
                                     mean = matrix(0, nrow = 3, ncol = n_eval),
                                     U = U, V=solve(ar1_cov(n_eval,ar_slope,1))
          )
          vec_1[obs,] <- realization[1,]
          vec_2[obs,] <- realization[2,]
          vec_3[obs,] <- realization[3,]
        }
        
        data_fun_1[[rep_number]] <- fda::Data2fd(argvals = eval_point,
                                                 y = t(vec_1),
                                                 basisobj = basis_predictor)
        data_fun_2[[rep_number]] <- fda::Data2fd(argvals = eval_point,
                                                 y = t(vec_2),
                                                 basisobj = basis_predictor)
        scalar.predictor <- matrix(NA, nrow = n_sample, ncol = n_predictor_scalar)
        for (ii in 1:n_sample){
          signal.mixture.chunk <-
            scalar.predictor[ii,] <- as.numeric(
              lapply(
                chunk2(
                  vec_3[ii,],
                  n_predictor_scalar
                ),
                mean)
            )
        }
        data_scalar[[rep_number]] <- scalar.predictor
        
        signal_functional[[rep_number]] <- inprod(reg_coef_1, data_fun_1[[rep_number]]) + inprod(reg_coef_2, data_fun_2[[rep_number]])
        signal_scalar[[rep_number]]  <- t(data_scalar[[rep_number]] %*% scalar_reg_coef)
      } # less than 10 seconds on laptop
      
      
      signal_multiplier <- (var(unlist(signal_functional)) / var(unlist(signal_scalar)))/9
      
      for (rep_number in 1:n_rep){
        data_scalar[[rep_number]]   <- sqrt(signal_multiplier) * data_scalar[[rep_number]]
        signal_scalar[[rep_number]] <- sqrt(signal_multiplier) * signal_scalar[[rep_number]]
        signal_total[[rep_number]] <- signal_functional[[rep_number]] +  signal_scalar[[rep_number]]
      }
      
      
      # generate response
      noise_var <- 1/9 * var(unlist(signal_total))
      for (rep_number in 1:n_rep){
        response[[rep_number]] <- signal_total[[rep_number]]+ sqrt(noise_var)*rnorm(n_sample, 0, 1)
      }
      
      save(response, data_fun_1, data_fun_2,data_scalar, response, response, eval_point, n_sample, n_basis_predictor, n_rep, lambda_pls_base, noise_var,file = "simul_data.RData")
      save.image(paste(datapath, paste0(trial_name, ".RData"), sep = "/")) # creating ".RData" in current working directory
      
      
      ```
      