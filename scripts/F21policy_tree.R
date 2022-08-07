### F21policy_tree.R
### Functions to fit policy tree algorithm and calculate a cross-validated value estimate of the resulting rule
### Separate functions for calculating CV value based on inner folds (for parameter tuning) vs. outer folds (evaluation of final rule)
### Working Directory = "WISDM-Precision-Medicine/"
### Author: Kushal Shah

library(policytree)
library(DiagrammeR)
library(dplyr)
select <- dplyr::select

#' Policy tree fitting and cross validated value estimate based on outer test folds
#' 
#' @param dat Sample size; needed if outer_set not provided
#' @param cv_folds Indices of outer folds corresponding to dataset (output of cv.folds() function)
#' @param K Number of outer folds
#' @param depth A scalar for the depth of the policy tree
#' 
#' @return Training and testing value estimates of fitted policy tree, along with values of CGM- and BGM-only rules, based on the specified outer folds
#' 
run.policytree.cv <- function(dat, cv_folds, K, depth) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val <- matrix(data = NA, nrow = K, ncol = 4)
  colnames(fold_val) <- c("train", "test", "cgm_test", "bgm_test")
  
  # Fitting algorithm on K outer training folds and estimating value on K outer testing folds
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
    bgm_value <- mean(dat_test_Y[dat_test_A == 0])
    
    fold_val[k, 1] <- train_value
    fold_val[k, 2] <- test_value
    fold_val[k, 3] <- cgm_value
    fold_val[k, 4] <- bgm_value
  }
  return(fold_val)
}

#' Policy tree fitting and cross validated value estimate based on inner test folds
#' 
#' @param dat Sample size; needed if outer_set not provided
#' @param cv_folds_inner Indices of inner folds corresponding to dataset
#' @param K Number of outer folds
#' @param L Number of inner folds
#' @param depth A scalar for the depth of the policy tree
#' 
#' @return Training and testing value estimates of fitted policy tree, along with values of CGM-only rule, based on the specified inner folds
#' 
run.policytree.cv.inner <- function(dat, cv_folds_inner, K, L, depth) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val_inner <- matrix(data = NA, nrow = K*L, ncol = 3)
  colnames(fold_val_inner) <- c("train_train", "train_test", "cgm_train_test")
  m = 0
  
  # Fitting algorithm on K*L inner training folds and estimating value on K*L inner testing folds
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
    }
  }
  return(fold_val_inner)
}
