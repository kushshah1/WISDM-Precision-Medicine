### C21policy_tree.R
### Script running functions for policy tree fitting and value estimation
### Working Directory = "WISDM-Precision-Medicine/"
### Author: Kushal Shah

source("./scripts/F10CV_base.R")
source("./scripts/F21policy_tree.R")
dat_clean_XAY_full <- readRDS("./data/dat_clean_XAY_full.rds")
dat <- dat_clean_XAY_full

# Outer CV folds setup
set.seed(1)
n <- nrow(dat)
K = 5 # Outer folds
L = K # Inner folds within each outer fold
cv_folds <- cv.folds(n, K)
train_folds <- cv_folds$train
test_folds <- cv_folds$test

# Inner CV folds setup
cv_folds_inner <- list(K)
for (k in 1:K) {
  cv_folds_inner[[k]] <- cv.folds(n = NA, K, is.inner = TRUE, outer_set = train_folds[[k]])  
}

# Policy Tree K-Fold Inner CV for Parameter Tuning
set.seed(2)

depth = 1
results_depth1 <- run.policytree.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)
colMeans(results_depth1)

depth = 2
results_depth2 <- run.policytree.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)
colMeans(results_depth2)

depth = 3
results_depth3 <- run.policytree.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)
colMeans(results_depth3)

# Policy Tree K-Fold Outer CV for Evaluation of Optimal Model
set.seed(2)
depth = 1
results <- run.policytree.cv(dat = dat, cv_folds = cv_folds, K = K, depth = depth)
colMeans(results) # Outer set training and testing (held-out) value estimate of optimal rule, along with value estimate of CGM-only and BGM-only rules on the same held-out test set

# SD of optimal rule, CGM-only, and BGM-only
sd(results[,2]) / sqrt(5) # Optimal
sd(results[,3]) / sqrt(5) # CGM-Only
sd(results[,4]) / sqrt(5) # BGM-Only

### Outputted final decision rule, trained on all data ###

  # 1. Preprocessing
  depth = 1
  nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
  dat_train_X_with_intercept <- model.matrix(~., dat %>% select(-all_of(nonFeatureVars)))
  dat_train_X <- dat_train_X_with_intercept[, 2:ncol(dat_train_X_with_intercept)]
  dat_train_Y <- as.matrix(dat$gluBelow70Chg)
  dat_train_A <- as.matrix(dat$TrtGroup)
  dat_train_A_replaceCGM <- replace(dat_train_A, dat_train_A == "CGM", 1)
  dat_train_A_replaceBGM <- replace(dat_train_A_replaceCGM, dat_train_A_replaceCGM == "BGM", 0)
  dat_train_A_numeric <- as.numeric(dat_train_A_replaceBGM)
  
  # 2. Visualization of Final Rule
  cf <- grf::causal_forest(X = dat_train_X, Y = dat_train_Y, W = dat_train_A_numeric)
  dr <- double_robust_scores(cf)
  tree.out = policy_tree(dat_train_X, dr, depth = depth)
  plot(tree.out) # Final Rule
  table(predict(tree.out, newdata = dat_train_X)) # Number of patients in each treatment group
  
  # 3. Outer Training Set Value of Optimal Rule
  dat_train_A_pred <- predict(tree.out, dat_train_X)
  dat_train_A_pred <- replace(dat_train_A_pred, dat_train_A_pred == 2, "CGM")
  dat_train_A_pred <- replace(dat_train_A_pred, dat_train_A_pred == 1, "BGM")
  train_value <- mean(dat_train_Y[dat_train_A == dat_train_A_pred])


# Creating dataset with column of optimal predicted treatments added
dat_with_opt_trt <- dat
dat_with_opt_trt$opt <- predict(tree.out, newdata = dat_train_X)
#saveRDS(dat_with_opt_trt, "./data/dat_with_opt_trt.rds")