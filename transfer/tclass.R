#####
## TRANSFER CLASSIFICATION
#####


library(magrittr)
library(Rcpp)
library(parallel)
library(glmnet)
library(class)

sourceCpp('ident_regression/data/pdist.cpp') # code from http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/
ddir <- "~/stat312data"
list.files(ddir)


load(paste0(ddir, "/voxels_train.RData"))
load(paste0(ddir, "/indexTrain.RData"))
index <- train_index
lookup <- matrix(0, 1750, 2)
for (i in 1:1750) lookup[i, ] <- which(index == i)
#whichrep <- index
#whichrep[lookup[, 1]] <- 1
#whichrep[lookup[, 2]] <- 2
#rownames(voxels) <- paste0("s", index, "_", whichrep)
#rownames(voxels)[1:10 + 2000]
#colnames(voxels) <- paste0("v", 1:25927)
#save(voxels, file = paste0(ddir, "/voxels_train.RData"))

nafilt <- voxels %>% colSums %>% is.na %>% `!`
voxels_c <- voxels[, nafilt]
nvox0 <- dim(voxels_c)[2]

## SNR analysis

withinss <- (voxels_c[lookup[, 1], ] - voxels_c[lookup[, 2], ])^2 %>% colSums
vars <- apply(voxels_c, 2, var)
res <- hist(withinss, breaks = 40)
histmat <- cbind(res$mids, res$counts)
# truncated hist mat
inds_trunc <- cumsum(histmat[, 2]) %>% { which( (. > nvox0 * .02) & (. < nvox * .98) ) } 
t_histmat <- histmat[inds_trunc, ]
plot(histmat[, 1], log(histmat[, 2]))
points(t_histmat[, 1], log(t_histmat[, 2]), pch = "+")
t_histmat <- cbind(t_histmat[, 1], t_histmat[, 1]^2, t_histmat[, 2])
res_pr <- glm(t_histmat[, 3] ~ t_histmat[, 1] + t_histmat[, 2], family = "poisson")
cfs <- res_pr$coefficients
pr_dens <- (cfs[1] + cfs[2] * histmat[, 1] + cfs[3] * histmat[, 1]^2) %>% exp %>% {. /sum(.)}
nrejs <- cumsum(histmat[, 2])
erejs <- cumsum(nvox0 * pr_dens)
fdp_est <- erejs/nrejs
plot(histmat[, 1], fdp_est, main = "estimated FDP", xlab = "withinss")
thres <- fdp_est %>% `<`(.2) %>% which %>% max %>% { histmat[., 1] }
thres # 2775
snr_filt <- withinss < thres

voxels_s <- voxels_c[, snr_filt]
(nvox <- dim(voxels_s)[2]) # 478

ivoxels <- cbind(index, voxels_s)

## PERFORMANCE ON A TEST SET GIVEN AN ESTIMATE OF 'A'

nclass <- 10
a <- diag(rep(1, nvox))
err_test <- function(te_set1, te_set2, a, nclass, nreps = 10) {
  n_te <- dim(te_set1)[1]
  res <- numeric(nreps)
  for (i in 1:nreps) {
    ch_inds <- sample(n_te, nclass)
    mu_a <- te_set1[ch_inds, ] %*% a
    y_a <- te_set2[ch_inds, ] %*% a
  }
}

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



