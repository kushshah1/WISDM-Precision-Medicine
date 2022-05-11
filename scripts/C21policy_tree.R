source("./scripts/F10CV_base.R")
source("./scripts/F21policy_tree.R")
dat_clean_XAY_full <- readRDS("./data/dat_clean_XAY_full.rds")

# CV folds setup
set.seed(1)
n <- nrow(dat_clean_XAY_full)
K = 5 # outer folds
L = K
cv_folds <- cv.folds(n, K)
train_folds <- cv_folds$train
test_folds <- cv_folds$test

# Inner CV folds setup
cv_folds_inner <- list(K)
for (k in 1:K) {
  cv_folds_inner[[k]] <- cv.folds(n = NA, K, is.inner = TRUE, outer_set = train_folds[[k]])  
}

# Checks
# c(cv_folds_inner[[1]]$test[[1]], cv_folds_inner[[1]]$train[[1]]) %in% train_folds[[1]] # all true
# c(cv_folds_inner[[1]]$test[[2]], cv_folds_inner[[1]]$train[[2]]) %in% train_folds[[1]] # all true
# c(cv_folds_inner[[2]]$test[[1]], cv_folds_inner[[2]]$train[[1]]) %in% train_folds[[1]] # shouldn't all be true

# Policy Tree K-Fold Inner CV for Parameter Tuning
set.seed(2)

depth = 1
results_depth1 <- run.policytree.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)
colMeans(results_depth1)

depth = 2
results_depth2 <- run.policytree.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)
colMeans(results_depth2)

depth = 3
results_depth3 <- run.policytree.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)
colMeans(results_depth3)

# Policy Tree K-Fold CV
set.seed(2)
depth = 2
results <- run.policytree.cv(dat = dat_clean_XAY_full, cv_folds = cv_folds, K = K, depth = depth)
colMeans(results)
mean((dat_clean_XAY_full %>% filter(TrtGroup == "CGM"))$gluBelow70Chg)






# Trained on all data
depth = 2
nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
dat_train_X_with_intercept <- model.matrix(~., dat_clean_XAY_full %>% select(-all_of(nonFeatureVars)))
dat_train_X <- dat_train_X_with_intercept[, 2:ncol(dat_train_X_with_intercept)]
dat_train_Y <- as.matrix(dat_clean_XAY_full$gluBelow70Chg)
dat_train_A <- as.matrix(dat_clean_XAY_full$TrtGroup)
dat_train_A_replaceCGM <- replace(dat_train_A, dat_train_A == "CGM", 1)
dat_train_A_replaceBGM <- replace(dat_train_A_replaceCGM, dat_train_A_replaceCGM == "BGM", 0)
dat_train_A_numeric <- as.numeric(dat_train_A_replaceBGM)

cf <- grf::causal_forest(X = dat_train_X, Y = dat_train_Y, W = dat_train_A_numeric)
dr <- double_robust_scores(cf)
tree.out = policy_tree(dat_train_X, dr, depth = depth)
plot(tree.out)
table(predict(tree.out, newdata = dat_train_X))

dat_train_A_pred <- predict(tree.out, dat_train_X)
dat_train_A_pred <- replace(dat_train_A_pred, dat_train_A_pred == 2, "CGM")
dat_train_A_pred <- replace(dat_train_A_pred, dat_train_A_pred == 1, "BGM")

train_value <- mean(dat_train_Y[dat_train_A == dat_train_A_pred])