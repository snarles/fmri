#####
## TRANSFER CLASSIFICATION
#####


library(magrittr)
library(Rcpp)
library(parallel)
library(glmnet)

sourceCpp('ident_regression/data/pdist.cpp') # code from http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
ddir <- "~/stat312data"
list.files(ddir)

nafilt <- voxels %>% colSums %>% is.na %>% `!`
voxels_c <- voxels[, nafilt]

load(paste0(ddir, "/voxels_train.RData"))
load(paste0(ddir, "/indexTrain.RData"))
index <- train_index
lookup <- matrix(0, 1750, 2)
for (i in 1:1750) lookup[i, ] <- which(index == i)

ivoxels <- cbind(index, voxels_c)

## GENERATE TRAINING SET, ETC

n_tr <- 1000
n_te <- 1750 - n_tr
tr_inds <- sample(1750, n_tr)
te_pick <- rbinom(n_te, 1, .5) + 1
te_inds <- setdiff(1:1750, tr_inds)

tr_set1 <- ivoxels[lookup[tr_inds, 1], ]
tr_set2 <- ivoxels[lookup[tr_inds, 2], ]
te_set1 <- ivoxels[lookup[cbind(te_inds, te_pick)], ]
te_set2 <- ivoxels[lookup[cbind(te_inds, 3 - te_pick)], ]


