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




X <- feature_train[converted_indices, ]
Yall <- train_resp_all[, filt1]
Yv1 <- Yall[, roi3$v1]
dim(Yv1) # 3500 1294
Yv2 <- Yall[, roi3$v2]
dim(Yv2) # 3500 2083
Yv3 <- Yall[, roi3$v3]
dim(Yv3) # 3500 1790
Yv4 <- Yall[, roi3$v4]
dim(Yv4) # 3500 1535


####
## Training and test partitions preserve pairs
####

Y <- Yv1[, sample(1294, 500)]
n_tr <- 200
n_te <- 100
s_ <- sample(1750, 1750)
s_tr <- s_[1:n_tr]
s_te <- s_[n_tr + (1:n_te)]
inds_tr <- which(converted_indices %in% s_tr)
inds_te <- which(converted_indices %in% s_te)
X_tr <- X[inds_tr, ]
X_te <- X[inds_te, ]
Y_tr <- Y[inds_tr, ]
Y_te <- Y[inds_te, ]

mcc <- 3
obs <- list(X = X_tr, Y = Y_tr, X_te = X_te)
p_CV <- do.call2(params_CV1, obs, filtr = FALSE, mc.cores = mcc, rule = "lambda.min")
pre_CV <- do.call(pre_mle, c(obs, p_CV))
CV_cl <- do.call(post_probs, c(obs, pre_CV))$cl
(CV_err <- sum(CV_cl != truth$i_chosen))


