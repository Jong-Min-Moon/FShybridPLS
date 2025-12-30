# Generated from create-FSHybridPLS.Rmd: do not edit by hand

#' Generate K-Fold Cross-Validation Indices
#'
#' Creates a list of training and validation indices for K-fold cross-validation.
#' Uses `caret::createFolds` to ensure stratified or random splitting.
#'
#' @param n_sample Integer. The total number of samples.
#' @param n_fold Integer. The number of folds.
#'
#' @return A list of lists, where each sub-list contains `idx_train` and `idx_valid`.
#' @export
create_idx_kfold <- function(n_sample, n_fold) {
  # Generate folds using caret (returns validation indices by default)
  idx_list <- caret::createFolds(seq_len(n_sample), k = n_fold, list = TRUE, returnTrain = FALSE)
  
  idx_list_list <- vector("list", length(idx_list))
  
  for (i in seq_along(idx_list)) {
    valid_hyper_index <- idx_list[[i]]
    # Training set is the complement of the validation set
    train_hyper_index <- setdiff(seq_len(n_sample), valid_hyper_index)
    
    idx_list_list[[i]] <- list(
      "idx_train" = train_hyper_index, 
      "idx_valid" = valid_hyper_index
    )
  }
  
  return(idx_list_list)
}
