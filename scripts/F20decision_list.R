### F20decision_list.R
### Functions to fit decision list algorithm and calculate a cross-validated value estimate of the resulting rule
### Separate functions for calculating CV value based on inner folds (for parameter tuning) vs. outer folds (evaluation of final rule)
### Working Directory = "WISDM-Precision-Medicine/"
### Author: Kushal Shah

library(listdtr)
library(dplyr)
select <- dplyr::select

#' Decision list fitting and cross validated value estimate based on outer test folds
#' 
#' @param dat Sample size; needed if outer_set not provided
#' @param cv_folds Indices of outer folds corresponding to dataset (output of cv.folds() function)
#' @param K Number of outer folds
#' @param maxlen A scalar for the maximum length of the decision list
#' 
#' @return Training and testing value estimates of fitted decision list, along with values of CGM- and BGM-only rules, based on the specified outer folds
#' 
run.decisionlist.cv <- function(dat, cv_folds, K, maxlen) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val <- matrix(data = NA, nrow = K, ncol = 4)
  colnames(fold_val) <- c("train", "test", "cgm_test", "bgm_test")
  
  # Fitting algorithm on K outer training folds and estimating value on K outer testing folds
  for (k in 1:K) {
    # 1. Preprocessing
    dat_train_Y <- as.matrix(dat$gluBelow70Chg[train_folds[[k]]])
    dat_train_A <- as.matrix(dat$TrtGroup[train_folds[[k]]])
    dat_train_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[train_folds[[k]], ]
    dat_train_X <- dat_train_X_with_intercept[, 2:ncol(dat_train_X_with_intercept)]
    
    dat_test_Y <- as.matrix(dat$gluBelow70Chg[test_folds[[k]]])
    dat_test_A <- as.matrix(dat$TrtGroup[test_folds[[k]]])
    dat_test_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[test_folds[[k]], ]
    dat_test_X <- dat_test_X_with_intercept[, 2:ncol(dat_test_X_with_intercept)]
    
    # 2. Modeling
    dtr.out <- listdtr(y = dat_train_Y, a = dat_train_A, x = dat_train_X, stage.x = rep(1, ncol(dat_train_X)), maxlen = maxlen)
    
    # 3. Evaluation
    dat_train_A_pred <- predict(dtr.out, dat_train_X, stage = 1)
    dat_test_A_pred <- predict(dtr.out, dat_test_X, stage = 1)
    
    train_value <- mean(dat_train_Y[dat_train_A == dat_train_A_pred])
    test_value <- mean(dat_test_Y[dat_test_A == dat_test_A_pred])
    cgm_value <- mean(dat_test_Y[dat_test_A == "CGM"])
    bgm_value <- mean(dat_test_Y[dat_test_A == "BGM"])
    
    print(k)
    fold_val[k, 1] <- train_value
    fold_val[k, 2] <- test_value
    fold_val[k, 3] <- cgm_value
    fold_val[k, 4] <- bgm_value
  }
  return(fold_val)
}

#' Decision list fitting and cross validated value estimate based on inner test folds
#' 
#' @param dat Sample size; needed if outer_set not provided
#' @param cv_folds_inner Indices of inner folds corresponding to dataset
#' @param K Number of outer folds
#' @param L Number of inner folds
#' @param maxlen A scalar for the maximum length of the decision list
#' 
#' @return Training and testing value estimates of fitted decision list, along with values of CGM-only rule, based on the specified inner folds
#' 
run.decisionlist.cv.inner <- function(dat, cv_folds_inner, K, L, maxlen) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val_inner <- matrix(data = NA, nrow = K*L, ncol = 3)
  colnames(fold_val_inner) <- c("train_train", "train_test", "cgm_train_test")
  m = 0
  
  # Fitting algorithm on K*L inner training folds and estimating value on K*L inner testing folds
  for (k in 1:K) {
    for (l in 1:L) {
      # 1. Preprocessing
      dat_train_Y <- as.matrix(dat$gluBelow70Chg[cv_folds_inner[[k]]$train[[l]]])
      dat_train_A <- as.matrix(dat$TrtGroup[cv_folds_inner[[k]]$train[[l]]])
      dat_train_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[cv_folds_inner[[k]]$train[[l]], ]
      dat_train_X <- dat_train_X_with_intercept[, 2:ncol(dat_train_X_with_intercept)]
      
      dat_test_Y <- as.matrix(dat$gluBelow70Chg[cv_folds_inner[[k]]$test[[l]]])
      dat_test_A <- as.matrix(dat$TrtGroup[cv_folds_inner[[k]]$test[[l]]])
      dat_test_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))[cv_folds_inner[[k]]$test[[l]], ]
      dat_test_X <- dat_test_X_with_intercept[, 2:ncol(dat_test_X_with_intercept)]
      
      # 2. Modeling
      dtr.out <- listdtr(y = dat_train_Y, a = dat_train_A, x = dat_train_X, stage.x = rep(1, ncol(dat_train_X)), maxlen = maxlen)
      
      # 3. Evaluation
      dat_train_A_pred <- predict(dtr.out, dat_train_X, stage = 1)
      dat_test_A_pred <- predict(dtr.out, dat_test_X, stage = 1)
      
      train_value <- mean(dat_train_Y[dat_train_A == dat_train_A_pred])
      test_value <- mean(dat_test_Y[dat_test_A == dat_test_A_pred])
      cgm_value <- mean(dat_test_Y[dat_test_A == "CGM"])
      
      m = m+1
      print(m)
      fold_val_inner[m, 1] <- train_value
      fold_val_inner[m, 2] <- test_value
      fold_val_inner[m, 3] <- cgm_value
    }
  }
  return(fold_val_inner)
}
