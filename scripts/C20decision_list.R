### C20decision_list.R
### Script running functions for decision list fitting and value estimation
### Working Directory = "WISDM-Precision-Medicine/"
### Author: Kushal Shah

source("./scripts/F10CV_base.R")
source("./scripts/F20decision_list.R")
dat_clean_XAY_full <- readRDS("./data/dat_clean_XAY_full.rds")
dat <- dat_clean_XAY_full

n <- nrow(dat)
K = 5 # Outer folds
L = K # Inner folds within each outer fold (required to match K in this script)

# Outer CV folds setup
set.seed(1)
cv_folds <- cv.folds(n, K)
train_folds <- cv_folds$train
test_folds <- cv_folds$test

# Inner CV folds setup
cv_folds_inner <- list(K)
for (k in 1:K) {
  cv_folds_inner[[k]] <- cv.folds(n = NA, K, is.inner = TRUE, outer_set = train_folds[[k]])  
}

# Decision List K-Fold Inner CV for Parameter Tuning
set.seed(2)
maxlen = 1L
results_len1 <- run.decisionlist.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len1)

set.seed(2)
maxlen = 2L
results_len2 <- run.decisionlist.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len2)

set.seed(2)
maxlen = 3L
results_len3 <- run.decisionlist.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len3)

set.seed(2)
maxlen = 4L
results_len4 <- run.decisionlist.cv.inner(dat = dat, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len4)