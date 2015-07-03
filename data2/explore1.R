####
## Load 3500 training data
####

library(magrittr)

ddir <- '/home/snarles/stat312data/'
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

converted_indices <- match(train_index, trainstim)
X <- feature_train[converted_indices, ]
filt_vox <- (!is.na(colSums(train_resp)))
filt_loc <- voxel.loc[filt_vox, ]
roi2 <- sapply(roi, function(v) { voxel.loc[v, ] }, USE.NAMES = TRUE)
roi2$v1 %>% head
Y <- t(train_resp[, filt_vox])
  
####
## Training and test partitions preserve pairs
####

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




