`%notin%` <- Negate(`%in%`)

cv.folds <- function(n = NA, K, is.inner = FALSE, outer_set = NA) {
  if (is.inner == FALSE) {
    row_indices <- 1:n
    data_random <- sample(row_indices)
  } else {
    data_random <- sample(outer_set)
    n <- length(outer_set)
  }
  
  rem = n %% K
  ceil <- ceiling(n/K)
  floor <- floor(n/K)
  check_n <- rem*ceil + (K-rem)*floor
  
  if (rem == 0) {
    fold_groupings <- rep((rem+1):K, each = floor)
  } else {
    fold_groupings <- c(rep(1:rem, each = ceil), rep((rem+1):K, each = floor)) 
  }
  
  test_folds_unsorted <- split(data_random, fold_groupings)
  test_folds <- lapply(test_folds_unsorted, sort, decreasing = FALSE)
  
  train_folds_unsorted <- list()
  for (i in 1:length(test_folds)) {
    train_folds_unsorted[[i]] <- data_random[data_random %notin% test_folds[[i]]]
  }
  train_folds <- lapply(train_folds_unsorted, sort, decreasing = FALSE)
  
  return(list(test = test_folds, train = train_folds))
}

#cv.folds.nested <- function(outer_folds) {
#  
#}

#cv_folds_2 <- cv.folds()
#outer_folds <- cv_folds
#outer <- outer_folds$train[[1]]
