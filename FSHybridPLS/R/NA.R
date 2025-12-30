# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' @title PFR Model Handler
#' @description R6 class to prepare data, fit, and evaluate a Penalized Functional Regression (PFR) model.
#'
#' @importFrom R6 R6Class
#' @importFrom refund pfr
#' @export
PFRModel <- R6::R6Class(
  "PFRModel",

  public = list(
    # --- Public Fields ---
    trainData = NULL,
    testData = NULL,
    model = NULL,
    n_scalar_predictor = NULL,
    n_basis = NULL,
    eval_point = NULL,
    
    # Prediction storage
    y_train_pred = NULL,
    y_test_pred = NULL,
    
    # Normalization factors for error reporting (optional usage)
    response_range_train = 1, # Default to 1 (standardized data usually implies scale ~1)
    response_range_test = 1,

    #' @title Initialize PFR Model
    #' @description Constructs the instance, splits/normalizes data, and formats it for `refund::pfr`.
    #'
    #' @param W A `predictor_hybrid` object containing all predictors.
    #' @param y Vector. Response variable.
    #' @param n_eval Integer. Number of evaluation points for functional predictors.
    #' @param train_ratio Numeric. Split ratio for training data.
    initialize = function(W, y, n_eval, train_ratio) {
      
      # 1. Run Preprocessing Pipeline (Split & Normalize)
      # Assumes `split_and_normalize.all` is available in the environment
      processed <- split_and_normalize.all(W, y, train_ratio)
      
      # 2. Set Metadata
      self$n_scalar_predictor <- W$n_scalar
      # Assume all functional predictors share the same basis size
      self$n_basis <- W$n_basis_list[1]
      self$eval_point <- W$eval_point # Use the points from the object directly
      
      # If n_eval is different from W$eval_point length, you might need interpolation.
      # Here we assume W$eval_point is sufficient, but if you strictly need a new grid:
      # eval_min <- min(W$eval_point)
      # eval_max <- max(W$eval_point)
      # self$eval_point <- seq(eval_min, eval_max, length.out = n_eval)

      # 3. Prepare Training Data
      # PFR works best with lists or data frames. We use a list to allow matrix columns easily.
      self$trainData <- list()
      
      # Add Scalar Predictors (X1, X2, ...)
      # We extract columns from Z and name them explicitly
      for(i in 1:self$n_scalar_predictor) {
        self$trainData[[paste0("X", i)]] <- processed$predictor_train$Z[, i]
      }
      
      self$trainData$response <- processed$response_train
      
      # Add Functional Predictors (F1, F2) converted to dense matrices (Sample x Time)
      # Transpose is needed because eval.fd returns (Time x Sample)
      self$trainData$F1 <- t(fda::eval.fd(self$eval_point, processed$predictor_train$functional_list[[1]]))
      self$trainData$F2 <- t(fda::eval.fd(self$eval_point, processed$predictor_train$functional_list[[2]]))
      
      # 4. Prepare Test Data
      self$testData <- list()
      for(i in 1:self$n_scalar_predictor) {
        self$testData[[paste0("X", i)]] <- processed$predictor_test$Z[, i]
      }
      self$testData$response <- processed$response_test
      self$testData$F1 <- t(fda::eval.fd(self$eval_point, processed$predictor_test$functional_list[[1]]))
      self$testData$F2 <- t(fda::eval.fd(self$eval_point, processed$predictor_test$functional_list[[2]]))
      
      # Store response ranges if needed for normalized error calculation
      # Since data is already standardized, ranges might be small.
      self$response_range_train <- diff(range(self$trainData$response))
      self$response_range_test  <- diff(range(self$testData$response))
    },

    #' @title Generate PFR Model Formula
    #' @description Constructs the formula string dynamically based on the selected option.
    #' Uses `presmooth = 'fpca.sc'` for functional terms.
    #'
    #' @param option Character. "all", "functional_only", or "scalar_only".
    #' @return A character string representing the PFR call.
    get_model_string = function(option = "all") {
      
      # Validate option
      if (!option %in% c("all", "functional_only", "scalar_only")) {
        stop("Error: 'option' must be one of 'all', 'functional_only', or 'scalar_only'.")
      }

      # 1. Build Scalar Predictor String
      scalar_part <- ""
      if (option != "functional_only") {
        # "X1 + X2 + ... + Xp"
        scalar_part <- paste0(paste0("X", 1:self$n_scalar_predictor), collapse = " + ")
        # Add leading " + " if needed later
        if (option == "all") scalar_part <- paste0(" + ", scalar_part)
      }

      # 2. Build Functional Predictor String
      functional_part <- ""
      if (option != "scalar_only") {
        # Note: self$eval_point is injected directly via string interpolation requires it to be available in scope
        # Or we rely on the data object. 'lf' usually needs 'argvals' passed explicitly if not in data.
        # Since we use eval(parse()), we reference 'self$eval_point'.
        
        term_f1 <- paste0(
          "lf(F1, argvals = self$eval_point, presmooth = 'fpca.sc', presmooth.opts = list(nbasis = ", 
          self$n_basis, "))"
        )
        term_f2 <- paste0(
          "lf(F2, argvals = self$eval_point, presmooth = 'fpca.sc', presmooth.opts = list(nbasis = ", 
          self$n_basis, "))"
        )
        functional_part <- paste(term_f1, term_f2, sep = " + ")
      }

      # 3. Combine into Full Formula
      # Result: "pfr(response ~ lf(...) + lf(...) + X1 + ..., data = self$trainData)"
      full_formula <- paste0("refund::pfr(response ~ ", functional_part, scalar_part, ", data = self$trainData)")
      
      return(full_formula)
    },

    #' @title Fit PFR Model
    #' @description Evaluates the formula string to fit the model.
    #' @param option Character. Model configuration option.
    fit = function(option = "all") {
      # Get the command string
      cmd <- self$get_model_string(option)
      
      # Execute the command in the current environment
      self$model <- eval(parse(text = cmd))
    },

    #' @title Compute Prediction Error
    #' @description Predicts on training and testing sets and computes normalized RMSE.
    #' @return A list with `train_error` and `test_error`.
    computeError = function() {
      if (is.null(self$model)) stop("Model not fitted. Run fit() first.")

      # Predict
      self$y_train_pred <- predict(self$model, newdata = self$trainData)
      self$y_test_pred  <- predict(self$model, newdata = self$testData)

      # Compute RMSE (Normalized by range, as per original code)
      # Note: Ensure ranges are non-zero to avoid division by zero
      r_train <- if(self$response_range_train > 0) self$response_range_train else 1
      r_test  <- if(self$response_range_test > 0)  self$response_range_test  else 1

      train_rmse <- sqrt(mean((self$trainData$response - self$y_train_pred)^2)) / r_train
      test_rmse  <- sqrt(mean((self$testData$response - self$y_test_pred)^2))  / r_test

      return(list(train_error = train_rmse, test_error = test_rmse))
    }
  )
)
