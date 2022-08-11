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

depth = 2
results_depth2 <- run.policytree.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)

depth = 3
results_depth3 <- run.policytree.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, depth = depth)

policy_tree_CV <- as.data.frame(rbind(c("Policy Tree", 1, colMeans(results_depth1)), c("Policy Tree", 2, colMeans(results_depth2)), 
                        c("Policy Tree", 3, colMeans(results_depth3))))
colnames(policy_tree_CV) <- c("Model", "Depth", "Training Set Value", "Inner Validation Set Value")
#saveRDS(policy_tree_CV, "./data/policy_tree_CV.rds")

# Policy Tree K-Fold Outer CV for Evaluation of Optimal Model
set.seed(2)
depth = 1
results <- run.policytree.cv(dat = dat, cv_folds = cv_folds, K = K, depth = depth)

values <- colMeans(results) # Outer set training/testing (held-out) value estimates of opt rule + value estimate of CGM-only and BGM-only rules on same held-out test set
SEs <- sapply(as.data.frame(results), sd) / sqrt(K) # SE of optimal rule (train/test), CGM-only, and BGM-only
optimal_CV <- rbind(values, SEs)
colnames(optimal_CV) <- c("Optimal Method - Training", "Optimal Method - Test",
                             "CGM-Only - Test", "BGM-Only - Test")
rownames(optimal_CV) <- c("Value Estimate", "SE")
#saveRDS(optimal_CV, "./data/optimal_CV.rds")


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
  
  # 2. Visualization of Final Rule (DiagrammeR used here)
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
dat_opt <- dat
dat_opt$opt <- predict(tree.out, newdata = dat_train_X)
#saveRDS(dat_opt, "./data/dat_opt.rds")