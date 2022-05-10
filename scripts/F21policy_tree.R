library(policytree)
library(DiagrammeR)
library(dplyr)
select <- dplyr::select

run.policytree.cv <- function(dat, cv_folds, K, depth) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val <- matrix(data = NA, nrow = K, ncol = 3)
  colnames(fold_val) <- c("train", "test", "cgm_test")
  for (k in 1:K) {
    # 1. Preprocessing
    dat_train_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[train_folds[[k]], ]
    dat_train_X <- dat_train_X_with_intercept[, 2:ncol(dat_train_X_with_intercept)]
    
    dat_train_A <- as.matrix(dat$TrtGroup[train_folds[[k]]])
    dat_train_A <- replace(dat_train_A, dat_train_A == "CGM", 1)
    dat_train_A <- replace(dat_train_A, dat_train_A == "BGM", 0)
    dat_train_A <- as.numeric(dat_train_A)
    
    dat_train_Y <- as.matrix(dat$gluBelow70Chg[train_folds[[k]]])
    
    dat_test_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[test_folds[[k]], ]
    dat_test_X <- dat_test_X_with_intercept[, 2:ncol(dat_test_X_with_intercept)]
    
    dat_test_A <- as.matrix(dat$TrtGroup[test_folds[[k]]])
    dat_test_A <- replace(dat_test_A, dat_test_A == "CGM", 1)
    dat_test_A <- replace(dat_test_A, dat_test_A == "BGM", 0)
    dat_test_A <- as.numeric(dat_test_A)
    
    dat_test_Y <- as.matrix(dat$gluBelow70Chg[test_folds[[k]]])
    
    # 2. Modeling
    cf <- grf::causal_forest(X = dat_train_X, Y = dat_train_Y, W = dat_train_A)
    dr <- double_robust_scores(cf)
    tree.out = policy_tree(dat_train_X, dr, depth = depth)
    
    # 3. Evaluation
    dat_train_A_pred <- predict(tree.out, dat_train_X) - 1 # Since predictions come as (1,2) but inputs were (0,1)
    dat_test_A_pred <- predict(tree.out, dat_test_X) - 1 # Since predictions come as (1,2) but inputs were (0,1)
    
    train_value <- mean(dat_train_Y[dat_train_A == dat_train_A_pred])
    test_value <- mean(dat_test_Y[dat_test_A == dat_test_A_pred])
    cgm_value <- mean(dat_test_Y[dat_test_A == 1])
    
    fold_val[k, 1] <- train_value
    fold_val[k, 2] <- test_value
    fold_val[k, 3] <- cgm_value
    print(fold_val)
  }
  return(fold_val)
}

run.policytree.cv.inner <- function(dat, cv_folds_inner, K, L, depth) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val_inner <- matrix(data = NA, nrow = K*L, ncol = 3)
  colnames(fold_val_inner) <- c("train_train", "train_test", "cgm_train_test")
  m = 0
  for (k in 1:K) {
    for (l in 1:L) {
      # 1. Preprocessing
      dat_train_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[cv_folds_inner[[k]]$train[[l]], ]
      dat_train_X <- dat_train_X_with_intercept[, 2:ncol(dat_train_X_with_intercept)]
      
      dat_train_A <- as.matrix(dat$TrtGroup[cv_folds_inner[[k]]$train[[l]]])
      dat_train_A <- replace(dat_train_A, dat_train_A == "CGM", 1)
      dat_train_A <- replace(dat_train_A, dat_train_A == "BGM", 0)
      dat_train_A <- as.numeric(dat_train_A)
      
      dat_train_Y <- as.matrix(dat$gluBelow70Chg[cv_folds_inner[[k]]$train[[l]]])
      
      dat_test_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[cv_folds_inner[[k]]$test[[l]], ]
      dat_test_X <- dat_test_X_with_intercept[, 2:ncol(dat_test_X_with_intercept)]
      
      dat_test_A <- as.matrix(dat$TrtGroup[cv_folds_inner[[k]]$test[[l]]])
      dat_test_A <- replace(dat_test_A, dat_test_A == "CGM", 1)
      dat_test_A <- replace(dat_test_A, dat_test_A == "BGM", 0)
      dat_test_A <- as.numeric(dat_test_A)
      
      dat_test_Y <- as.matrix(dat$gluBelow70Chg[cv_folds_inner[[k]]$test[[l]]])
      
      # 2. Modeling
      cf <- grf::causal_forest(X = dat_train_X, Y = dat_train_Y, W = dat_train_A)
      dr <- double_robust_scores(cf)
      tree.out = policy_tree(dat_train_X, dr, depth = depth)
      
      # 3. Evaluation
      dat_train_A_pred <- predict(tree.out, dat_train_X) - 1 # Since predictions come as (1,2) but inputs were (0,1)
      dat_test_A_pred <- predict(tree.out, dat_test_X) - 1 # Since predictions come as (1,2) but inputs were (0,1)
      
      train_value <- mean(dat_train_Y[dat_train_A == dat_train_A_pred])
      test_value <- mean(dat_test_Y[dat_test_A == dat_test_A_pred])
      cgm_value <- mean(dat_test_Y[dat_test_A == 1])
      
      m = m+1
      print(m)
      fold_val_inner[m, 1] <- train_value
      fold_val_inner[m, 2] <- test_value
      fold_val_inner[m, 3] <- cgm_value
      #print(fold_val_inner)
    }
  }
  return(fold_val_inner)
}
