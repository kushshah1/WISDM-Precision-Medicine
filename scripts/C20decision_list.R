source("./scripts/F10CV_base.R")
source("./scripts/F20decision_list.R")
dat_clean_XAY_full <- readRDS("./data/dat_clean_XAY_full.rds")

# CV folds setup
set.seed(1)
n <- nrow(dat_clean_XAY_full)
K = 5
L = K
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
results_len1 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len1)

set.seed(2)
maxlen = 2L
results_len2 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len2)

set.seed(2)
maxlen = 3L
results_len3 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len3)

set.seed(2)
maxlen = 4L
results_len4 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len4)

set.seed(2)
maxlen = 6L
results_len6 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len6)

set.seed(2)
maxlen = 8L
results_len8 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len8)

set.seed(2)
maxlen = 10L
results_len10 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len10)

set.seed(2)
maxlen = 12L
results_len12 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len12)

set.seed(2)
maxlen = 14L
results_len14 <- run.decisionlist.cv.inner(dat = dat_clean_XAY_full, cv_folds_inner = cv_folds_inner, K = K, L = L, maxlen = maxlen)
colMeans(results_len14)

# Decision List K-Fold CV
set.seed(2)
maxlen = 6L
results <- run.decisionlist.cv(dat = dat_clean_XAY_full, cv_folds = cv_folds, K = K, maxlen = maxlen)
colMeans(results)
mean((dat_clean_XAY_full %>% filter(TrtGroup == "CGM"))$gluBelow70Chg)


# Trained on all data
nonFeatureVars <- c("TrtGroup", "gluBelow70Chg")
train_Y <- as.matrix(dat_clean_XAY_full$gluBelow70Chg)
train_A <- as.matrix(dat_clean_XAY_full$TrtGroup)

train_X <- model.matrix(~., dat_clean_XAY_full %>% select(-all_of(nonFeatureVars)))
train_X <- train_X[, 2:ncol(train_X)]

set.seed(2)
dtr.out <- listdtr(train_Y, train_A, train_X, stage.x = rep(1, ncol(train_X)), maxlen = 6L)
print(dtr.out)
b <- table(predict(dtr.out, train_X, stage = 1))