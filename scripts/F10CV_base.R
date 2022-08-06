### F10CV_base.R
### Basic function to set up inner and outer training and testing folds for nested K-fold cross validation
### Working Directory = "WISDM-Precision-Medicine/"
### Author: Kushal Shah

`%notin%` <- Negate(`%in%`)

#' Indices for Nested K-Fold Cross Validation
#' 
#' @param n Sample size; needed if outer_set not provided
#' @param K Number of folds
#' @param is.inner Logical indicating whether inner folds are being created
#' @param outer_set If inner folds are to be created, the outer set must be provided as integer vector of indices
#' 
#' @return List of K integer vectors, each vector representing indices for one of K folds
#' 
cv.folds <- function(n = NA, K, is.inner = FALSE, outer_set = NA) {
  # If outer folds being created, create one index for each observation (1:n) and randomize the order
  if (is.inner == FALSE) {
    row_indices <- 1:n
    data_random <- sample(row_indices)
  # If inner folds being created, use the provided outer fold indices and randomize the order
  } else {
    data_random <- sample(outer_set)
    n <- length(outer_set)
  }
  
  # Internal check for fold sizes adding up to sample size of input data
  rem = n %% K
  ceil <- ceiling(n/K)
  floor <- floor(n/K)
  check_n <- rem*ceil + (K-rem)*floor
  if(n != check_n)
    warning("Warning: Sizes of CV folds do not add up to inputted sample size.")
  
  # Creates vector of length n indicating fold membership for each observation
  # Ensures floor(n/K) or ceiling(n/K) indices in each fold
  if (rem == 0) {
    fold_groupings <- rep((rem+1):K, each = floor)
  } else {
    fold_groupings <- c(rep(1:rem, each = ceil), rep((rem+1):K, each = floor)) 
  }
  
  # Splits indices into test folds according to fold membership specified in `fold_groupings`
  test_folds_unsorted <- split(data_random, fold_groupings)
  test_folds <- lapply(test_folds_unsorted, sort, decreasing = FALSE)
  
  # Creates training folds based on test folds
  train_folds_unsorted <- list()
  for (i in 1:length(test_folds)) {
    train_folds_unsorted[[i]] <- data_random[data_random %notin% test_folds[[i]]]
  }
  train_folds <- lapply(train_folds_unsorted, sort, decreasing = FALSE)
  
  return(list(test = test_folds, train = train_folds))
}
