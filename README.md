# FSHybridPLS: Functional and Scalar Hybrid Partial Least Squares

**FSHybridPLS** is an R package designed to perform Partial Least Squares (PLS) regression on "hybrid" predictors. A hybrid predictor is a single mathematical object containing both **functional data** (curves, time-series represented as `fd` objects) and **scalar covariates** (standard numeric matrices).

This package defines a joint Hilbert space $\mathcal{H} = \mathcal{F} \times \mathbb{R}^p$ and implements the arithmetic and algorithms necessary to perform penalized PLS directly on this space.
This package provides an R implementation of my paper, “Hybrid Partial Least Squares Regression with Multiple Functional and Scalar Predictors,” which will be released on arXiv soon. The work was co-authored with Professor Jeong Hoon Jang of the University of Texas Medical Branch.

## Installation

This package is generated using `litr`, a literate programming tool developed by Jacob Bien and Patrick Vossler (https://jacobbien.github.io/litr-project/). I used litr to build the package and you don't have to install litr just to use the package. To install the package:

```r
devtools::install("your_downloaded_directory/FSHybridPLS")
```



-----

## 1\. Core Data Structure: `predictor_hybrid`

The package revolves around the `predictor_hybrid` S3 class. This class acts as a container allowing mathematical operations to be performed on mixed data types simultaneously.

### Class Constructors

  * **`predictor_hybrid(Z, functional_list, eval_point)`**

      * **Description:** The main constructor. It calculates and stores the Gram matrices ($J$) and second-derivative penalty matrices ($\ddot{J}$) upon initialization for efficiency.
      * **Input:** \* `Z`: Numeric matrix (Scalars).
          * `functional_list`: List of `fd` objects (Functional).
          * `eval_point`: Vector of evaluation points.
      * **Output:** An object of class `predictor_hybrid`.

  * **`predictor_hybrid_from_coef(format, coef)`**

      * **Description:** Reconstructs a hybrid object from a flattened coefficient vector (used during optimization steps).
      * **Input:** A template `format` object and a numeric `coef` vector.

### Class Manipulation

  * **`subset_predictor_hybrid(W, i)`**
      * **Description:** Extracts the $i$-th observation (sample) as a new hybrid object.
  * **`replace_obs_hybrid(W, i, new_W)`**
      * **Description:** Replaces the $i$-th observation in `W` with `new_W`.

-----

## 2\. Hilbert Space Arithmetic

These functions implement the algebra required for the PLS algorithm. They handle the broadcasting of scalars and the integration of functional components.

  * **`inprod.predictor_hybrid(xi_1, xi_2)`**
      * **Description:** Computes the global inner product $\langle \mathbf{h}_1, \mathbf{h}_2\rangle_{\mathcal{H}}$. It sums the $L^2$ inner products of the functions and the dot product of the scalars.
  * **`inprod_pen.predictor_hybrid(xi_1, xi_2, lambda)`**
      * **Description:** Computes the penalized inner product using the second derivative of the functions, weighted by `lambda`.
  * **`add.predictor_hybrid(xi_1, xi_2, alpha)`**
      * **Description:** Computes $\xi_1 + \alpha \cdot \xi_2$. Supports broadcasting (e.g., adding a single mean hybrid object to a set of observations).
  * **`subtr.predictor_hybrid(input, other, alpha)`**
      * **Description:** Computes `input` - `alpha` \* `other`.
  * **`scalar_mul.predictor_hybrid(input, scalar)`**
      * **Description:** Scales all functional and scalar components by a numeric value.
  * **`add_broadcast(input, other)`** & **`subtr_broadcast(input, other)`**
      * **Description:** S4 methods to handle broadcasting arithmetic specifically for the matrix components.

-----

## 3\. Penalized PLS Algorithm (NIPALS)

These functions implement the iterative extraction of Partial Least Squares components.

### Main Solvers

  * **`nipals_pen_hybrid(W, y, n_iter, lambda, verbose)`**
      * **Description:** The primary user-facing algorithm. It runs the NIPALS loop: finding direction `xi`, calculating score `rho`, and residualizing `W` and `y`.
      * **Input:** Predictors `W`, response `y`, iterations `n_iter`, and smoothing parameters `lambda`.
  * **`fit.hybridPLS(W, y, n_iter, lambda)`**
      * **Description:** An internal wrapper for the iterative fitting process.
  * **`fit.hybridPLS_eigen(W, y, n_iter, lambda)`**
      * **Description:** A variation of the fitting algorithm using explicit eigen-decomposition of the constraint matrix.

### Component Computation

  * **`get_pls_comp(W, y, L)`**
      * **Description:** Solves the generalized eigenvalue problem $\max Cov(W, y)^2$ s.t. constraints.
      * **Input:** Residualized `W`, `y`, and Cholesky factor `L` of the constraint matrix.
      * **Output:** The PLS direction vector `xi_hat`.
  * **`get_xi_hat_linear_pen(W, y, lambda)`**
      * **Description:** Alternative formulation to solve for the PLS direction using a linear system solver.
  * **`get_rho(d_vec, W)`**
      * **Description:** Computes the PLS scores given the direction vector coefficients `d_vec`.

### Residualization Steps

These functions project the data onto the finding components and remove that information for the next iteration.

  * **`get_nu(y, rho)`**: Calculates the scalar regression coefficient for the response.
  * **`residualize_y(y, rho, nu)`**: Updates $y \leftarrow y - \nu \rho$.
  * **`get_delta(W, rho)`**: Calculates the hybrid regression coefficient (projection) for the predictors.
  * **`residualize_predictor(W, rho, delta)`**: Updates $W \leftarrow W - \rho \delta$.
  * **`LSE_hybrid`** & **`LSE_ptws`**: Least Squares Estimation helpers for the residualization steps.
  * **`get_pre_corrcov(W, y)`**: Computes the corrected cross-covariance matrix $V^*$.
  * **`get_CJ(W)`**: Helper to compute the product of coefficients and Gram matrices.

-----

## 4\. Matrix Construction (Regularization)

These functions build the specific block-diagonal matrices required to solve the penalized eigenvalue problem $J^* + \Lambda \ddot{J}^*$.

  * **`get_gram_matrix_block(obj)`**: Constructs the block-diagonal Gram matrix ($J^*$) of basis inner products.
  * **`get_penalty_hybrid(obj)`**: Constructs the penalty matrix ($\ddot{J}^*$) using second derivatives.
  * **`get_smoothing_param_hybrid(W, lambda)`**: Constructs the smoothing matrix $\Lambda$.
  * **`get_constraint_matrix_hybrid(W, lambda)`**: Combines the above into the final LHS matrix for the eigenproblem.

-----

## 5\. Simulation & Preprocessing

Tools for data generation, splitting, and normalization.

### Splitting

  * **`split.all(W_hybrid, response, train_ratio)`**: Splits hybrid data and response into train/test sets.
  * **`create_idx_train_test`** / **`get_idx_train`**: Helpers to generate random indices.
  * **`create_idx_kfold`**: Generates indices for K-fold cross-validation.

### Normalization

  * **`curve_normalize_train_test(train, test)`**: Centers and scales functional predictors. Crucially, it scales the *test* set using *training* statistics.
  * **`scalar_normalize_train_test(predictor_train, predictor_test)`**: Centers and scales scalar predictors.
  * **`btwn_normalize_train_test(train, test)`**: Balances the total variance between the functional block and the scalar block so neither dominates the optimization purely due to units.
  * **`split_and_normailize.all`**: A master wrapper that splits and applies all three normalizations.
  * **`curve_normalize`** / **`scalar_normalize`**: Low-level normalization logic.

### Helpers

  * **`is_same_basis(input, other)`**: Checks if two hybrid objects share the same functional basis.
  * **`compute_gram_matrix(basis)`**: Wrapper for `fda::inprod`.
  * **`are_all_gram_matrices_identical`**: Utility to ensure consistency.
  * **`rep_fd(fd_list, n)`**: Replicates an `fd` object `n` times (broadcasting).
  * **`n_sample.fd(fd_obj)`**: Returns the number of curves in an `fd` object.

-----

## 6\. Baseline: Modified MFPCA

The package includes a modified version of Multivariate Functional PCA (based on the `funData` package) to support simultaneous score prediction for training and test sets.

  * **`fit.mfpca_regression(W, y, ...)`**: Runs a baseline comparison model: MFPCA for dimension reduction $\to$ Linear Regression.
  * **`MFPCA`**: Modified main function allowing `mFData_predict` as input.
  * **`fpcaBasis`**: Modified wrapper for the PACE algorithm.
  * **`univDecomp`**: Dispatches univariate decompositions.
  * **`calcMFPCA`** / **`calcMFPCA_predict`**: Internal functions to calculate scores and eigenvectors.
  * **`calcBasisIntegrals`**: Computes basis integrals.
  * **`PACE`** / **`.PACE`**: Principal Analysis by Conditional Expectation implementation.
  * **`expandBasisFunction`** / **`univExpansion`**: Helpers to map scores back to function space.
