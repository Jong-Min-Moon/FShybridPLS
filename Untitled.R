




























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
