library(listdtr)
library(dplyr)
select <- dplyr::select

run.decisionlist.cv <- function(dat, cv_folds, K, maxlen) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val <- matrix(data = NA, nrow = K, ncol = 4)
  colnames(fold_val) <- c("train", "test", "cgm_test", "bgm_test")
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
    #print(fold_val)
  }
  return(fold_val)
}

run.decisionlist.cv.inner <- function(dat, cv_folds_inner, K, L, maxlen) {
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  fold_val_inner <- matrix(data = NA, nrow = K*L, ncol = 3)
  colnames(fold_val_inner) <- c("train_train", "train_test", "cgm_train_test")
  m = 0
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
      #print(fold_val_inner)
    }
  }
  return(fold_val_inner)
}
