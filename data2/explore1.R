####
## Load 3500 training data
####

library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
library(glmnet, warn.conflicts = FALSE)
library(prodlim)
source('utils/zattach.R')
source('transfer/source.R')
source('eb_ident/eigenprism.R')
source('eb_ident/bayes_reg.R')
source('eb_ident/source.R')


ddir <- '/home/snarles/stat312data/'
ddir <- '/home/rstudio/stat312data/'

list.files(ddir)
load(paste0(ddir, 'train_resp_all.RData'))
dim(train_resp) # 1750 25915
load(paste0(ddir, 'roi.RData'))
dim(voxel.loc) # 25915 3
feat_attr <- read.table(paste0(ddir, "featAttr.csv"), header = TRUE, sep = ",")
dim(feat_attr)
load(paste0(ddir, 'photos_stimindex.RData'))
train_index <- read.table(paste0(ddir, "indexTrain.csv"), header = FALSE, sep = ",") %>% as.numeric
load(paste0(ddir, 'feature_train.RData'))
dim(feature_train) # 1750 10921
#train_resp_all <- read.csv('~/stat312data/allVoxTrain.csv', header = FALSE, sep = ",", na.strings = "NaN")
train_resp_all <- t(read.csv(gzfile('~/stat312data/allVoxTrain.csv.gz'),
                           header = FALSE, sep = ",", na.strings = c("NA", "NaN", "N"),
                           stringsAsFactors = FALSE))
dim(train_resp_all) #   3500 25927


converted_indices <- match(train_index, trainstim)
i_set1 <- match(1:1750, converted_indices)
i_set2 <- 3501 - match(1:1750, rev(converted_indices))
length(unique(c(i_set1, i_set2)))

train_resp_1 <- train_resp_all[i_set1, ]
train_resp_2 <- train_resp_all[i_set2, ]

dim(train_resp_1) # 1750 25927

filt0 <- !is.na(colSums(train_resp))
filt1 <- !is.na(colSums(train_resp_1))
filt2 <- !is.na(colSums(train_resp_2))

train_resp_f <- train_resp[, filt0]
train_resp_1f <- train_resp_1[, filt1]
train_resp_2f <- train_resp_2[, filt2]

dim(train_resp_f) # 1750 22733
dim(train_resp_1f) # 1750 22733
dim(train_resp_2f) # 1750 22733

filt_loc <- voxel.loc[filt0, ]
roi2 <- sapply(roi, function(v) { voxel.loc[v, ] }, USE.NAMES = TRUE)
roi3 <- sapply(roi2, function(v) {
  row.match(as.data.frame(v), filt_loc) %>% {.[!is.na(.)]}
}, USE.NAMES = TRUE)

#samp <- sample(22733, 100)
#image(cor(train_resp_f[, samp], train_resp_2f[, samp]))




Xall <- feature_train[converted_indices, ]
Yall <- train_resp_all[, filt1]
Yv1 <- Yall[, roi3$v1]
dim(Yv1) # 3500 1294
Yv2 <- Yall[, roi3$v2]
dim(Yv2) # 3500 2083
Yv3 <- Yall[, roi3$v3]
dim(Yv3) # 3500 1790
Yv4 <- Yall[, roi3$v4]
dim(Yv4) # 3500 1535
X <- Xall[, 1:5000]
FF <- feature_train[, 1:5000]

####
## New functions
####

sscore <- function(ll, i_chosen) {
  mat <- ll$pprobs
  mat <- apply(mat, 1,  function(v) {
    v <- v - max(v)
    v <- exp(v)
    v/sum(v)
  })
  sum(t(mat)[cbind(1:dim(mat)[1], i_chosen)])
}
topk <- function(ll, i_chosen, k = 10) {
  mat <- ll$pprobs
  mmat <- cbind(i_chosen, mat)
  res <- apply(mmat, 1, function(v) v[1] %in% order(-v[-1])[1:k])
  sum(res)
}

####
## Training and test partitions preserve pairs
####

inds_Y <- sample(1294, 500)
Y <- Yv1[, inds_Y]
n_tr <- 200
n_te <- 100
s_ <- sample(1750, 1750)
s_tr <- sort(s_[1:n_tr])
s_te <- sort(s_[n_tr + (1:n_te)])
inds_tr <- which(converted_indices %in% s_tr)
inds_te <- which(converted_indices %in% s_te)
X_tr <- X[inds_tr, ]
X_te <- FF[s_te, ]
Y_tr <- Y[inds_tr, ]
Y_te <- Y[inds_te, ]
i_chosen <- match(converted_indices[inds_te], s_te)
pX <- dim(X)[2]

mcc <- 39
obs <- list(X = X_tr, Y = Y_tr, X_te = X_te, y_star = Y_te)
t1 <- proc.time()
p_CV <- do.call2(params_CV1, obs, 
                 filtr = FALSE, mc.cores = mcc, 
                 rule = "lambda.1se")
proc.time() - t1
t1 <- proc.time()
pre_CV <- do.call(pre_mle, c(obs, p_CV))
proc.time() - t1
t1 <- proc.time()
CV_cl <- do.call(post_probs, c(obs, pre_CV))
proc.time() - t1
(CV_good <- sscore(CV_cl, i_chosen))
(CV_crr <- sum(CV_cl$cl == i_chosen))
(CV_crr10 <- topk(CV_cl, i_chosen, 10))

B <- p_CV$B
T2 <- colSums(B^2) - 1e-5
t1 <- proc.time()
pre_EB <- do.call2(predictive_EB, c(obs, p_CV), T2 = T2, mc.cores = mcc)
proc.time() - t1
t1 <- proc.time()
EB_cl <- do.call(post_probs, c(obs, pre_EB))
proc.time() - t1

(EB_good <- sscore(EB_cl, i_chosen))
(EB_crr <- sum(EB_cl$cl == i_chosen))
(EB_crr10 <- topk(EB_cl, i_chosen, 10))
