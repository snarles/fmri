#############################################################
##                DATA PREPROCESSING                       ##
#############################################################

library(magrittr)
isqrtm <- function(m) {
  res <- eigen(m)
  d <- res$values
  if (min(d) < -1e-5) warning("Negative eigenvalues in isqrtm")
  d[d < 0] <- 0
  d[d > 0] <- 1/sqrt(d[d > 0])
  v <- res$vectors
  return (v %*% diag(d) %*% t(v))
}


library(Rcpp)
sourceCpp('ident_regression/data/pdist.cpp') # code from http://blog.felixriedel.com/2013/05/pairwise-distances-in-r/

ddir <- "~/stat312data"
#ddir <- "/home/ubuntu/stat312data"
list.files(ddir)

#load(paste0(ddir, "/all_voxel_locations.RData"))
#dim(voxel.loc) # 25915 3
load(paste0(ddir, "/roi.RData"))

#temp <- read.csv(paste0(ddir, "/allVoxTrain.csv"), header = FALSE,
#                 stringsAsFactors = FALSE)
#temp[, 1] <- as.numeric(temp[, 1])
#table(sapply(temp, class))
#voxels <- t(temp)
#dim(voxels)
dim(voxel.loc)
cnv <- character(length(train_index))
for (i in 1:length(roi)) {
  nn <- roi %>% names %>% `[[`(i)
  cnv[roi[[i]]] <- paste0(nn, "_", 1:length(roi[[i]]))
}
#save(voxels, file = paste0(ddir, "/voxels_train.RData"))

load(paste0(ddir, "/voxels_train.RData"))

train_index <- read.csv(paste0(ddir, "/indexTrain.csv"), header = FALSE)
train_index <- as.numeric(train_index)
u_inds <- unique(train_index)
length(u_inds)
temp <- numeric(max(u_inds))
temp[u_inds] <- 1:1750
train_index <- temp[train_index]
save(train_index, file = paste0(ddir, "/indexTrain.RData"))
load(paste0(ddir, "/indexTrain.RData"))

load(paste0(ddir, "/feature_train.RData"))

dim(feature_train) # 1750 10921
#feature_train2 <- feature_train[train_index, ]


#train_resp <- read.csv(paste0(ddir, "/train_resp_all.csv"), header = FALSE)
#dim(train_resp) #25915 1750
#train_resp <- t(train_resp)
#colnames(train_resp) <- colnames(voxels)
#save(train_resp, file = paste0(ddir, "/train_resp_all.RData"))
load(paste0(ddir, "/train_resp_all.RData"))
dim(train_resp)


## create the per-image average

train_avg <- matrix(0, 1750, dim(voxels)[2])
for (i in 1:1750) train_avg[i, ] <- voxels[train_index == i, ] %>% colMeans

nacount <- train_avg %>% apply(2, function(v) sum(is.na(v)))
nacount2 <- train_resp %>% apply(2, function(v) sum(is.na(v)))
table(nacount)
table(nacount2)
match_counts <- numeric()
for (i in unique(nacount2)) {
  if (sum(nacount2== i) == sum(nacount == i)) match_counts <- c(match_counts, i)
}
match_counts

nafilt <- nacount %in% match_counts
nafilt2 <- nacount2 %in% match_counts
sum(nafilt)
sum(nafilt2)


train_avg_filt <- train_avg[, nafilt]
train_resp_filt <- train_resp[, nafilt2]

i <- 2
voxels[train_index ==i, ] %>% t %>% plot
voxcc <- voxels %>% colSums %>% is.na %>% `!`
diffs <- sapply(1:1750, function(v) {
  temp <- voxels[train_index == i, voxcc]
  sum((temp[1, ] - temp[2, ])^2)
})
mean(diffs)
x <- sample(3500, 1750); (voxels[x, voxcc] - voxels[-x, voxcc]) %>% `^`(2) %>% rowSums %>% mean




cc <- train_avg_filt %>% colSums %>% is.na
train_avg_cc <- train_avg_filt[, !cc]
train_resp_cc <- train_resp_filt[, !cc]

avg_vars <- train_avg_cc %>% apply(1, var)
resp_vars <- train_resp_cc %>% apply(1, var)
plot(avg_vars, resp_vars)


dim(train_avg_cc)
i <- 20000
plot(train_avg_cc[, i], train_resp_cc[, i])
i <- 50
plot(train_avg_cc[i, ], train_resp_cc[i, ])

(train_avg_cc - train_resp_cc)^2 %>% rowSums %>% mean
(train_avg_cc - train_resp_cc[sample(1750, 1750)])^2 %>% rowSums %>% mean



## verify

for (i in 1:1750) {
  
}


## no obvious correlation effect

gaps <- sapply(1:1750, function(i) { which(train_index == i) %>% {max(.) - min(.)}  })
intradist <- sapply(1:1750, function(i) {
    temp <- voxels[train_index == i, complete_vox]
    sum((temp[, 1] - temp[, 2])^2)
  })
plot(gaps, intradist)
cfs <- lm(intradist ~ gaps)$coefficients
abline(cfs[1], cfs[2], col = "red")
mean(intradist)

plot(gaps, intradist, pch = ".")
points(1:132, sapply(1:132, function(i) mean(intradist[gaps == i])), cex = 3, col = "blue")

feat_attr <- read.csv(paste0(ddir, "/featAttr.csv"), header = TRUE,
                 stringsAsFactors = FALSE)
feat_lv <- feat_attr[2, ]

#####
## PROCESSING IMAGE FEATURES
#####

inds_train <- 1:1750
inds_valid <- 1750 + 1:120
features_all <- rbind(feature_train, feature_valid)
vars <- apply(features_all, 2, var)
lvars <- log(apply(features_all, 2, var))
plot(sort(lvars), type ="l")
var_filt <- (lvars > -10)
sum(var_filt)
dim(feat_attr)

comp_var <- sapply(
  1:4,
  function(i) {
    temp_filt <- var_filt & (feat_lv == i)
    median(vars[temp_filt])
  })
comp_var

for (i in 1:4) {
  features_all[, feat_lv == i] <-
    features_all[, feat_lv == i]/sqrt(comp_var[i])
}

features_all <- features_all[, var_filt]
features_train <- features_all[inds_train, ]
features_valid <- features_all[inds_valid, ]
feat_attr <- feat_attr[, var_filt]

dim(features_train)
train_index
#x_train <- features_train[train_index, ]
dim(train_v1)
length(train_index)
length(unique(train_index))
max(train_index)
max(valid_index)

####
## COVARIANCE OF ERROR
####

dim(features_train)
dim(train_v1)
train_index[1:10]
dm_train_v1 <- train_v1
i <- train_index[1]
for (i in unique(train_index)) {
    filt <- train_index == i
    dm_train_v1[, filt] <-
      t(apply(train_v1[, filt], 1, function(v) v - mean(v)))
}
dim(dm_train_v1)
sigma_e <-cov(t(dm_train_v1))
eye <- mean(diag(sigma_e)) * diag(rep(1, 100))
sigma_e <- 0.5 * sigma_e + 0.5 * eye
omega_e <- isqrtm(sigma_e)

####
## REGRESSION
####

library(parallel)
library(glmnet)
#cl <- makeCluster(5)

dim(features_train)
dim(train_resp)

lambdas <- 0:10/10
nlambdas <- length(lambdas)

prfunc <- function(i) {
    as.numeric(train_resp[1,])
    res <- glmnet(features_train, as.numeric(train_resp[i, ]), standardize = FALSE)
    pr <- predict(res, features_valid, s=lambdas)
    pr
}

res <- lapply(1:100, prfunc)

pr_error <- numeric(nlambdas)
misc_error <- numeric(nlambdas)
for (i in 1:nlambdas) {
  pvalid <- matrix(0, 120, 100)
  for (j in 1:100) {
    pvalid[, j] <- res[[j]][, i]
  }
  yhat <- pvalid[valid_index, ]
  ys <- t(valid_v1)
  pr_error[i] <- sum((yhat - ys)^2)
  for (z in 1:120) {
    y <- apply(valid_v1[, valid_index == z], 1, mean)
    diff <- t(pvalid) - y # 100 120
    cdiff <- omega_e %*% diff
    ds <- apply(cdiff^2, 2, sum)
    zhat <- order(ds)[1]
    misc_error[i] <- misc_error[i] + (zhat != z)
  }
}

plot(lambdas, misc_error)
plot(lambdas, pr_error)


####
## REGRESSION : USING TRAIN_RESP
####

lambdas <- 0:1000/40000
nlambdas <- length(lambdas)

ntrials <- 200
misc_errors <- matrix(0, ntrials, nlambdas)
pr_errors <- matrix(0, ntrials, nlambdas)

library(class)

proc.time()
for (ii in 1:ntrials) {
    tr_inds <- sample(1750, 1725, FALSE)
    te_inds <- setdiff(1:1750, tr_inds)
    nte <- length(te_inds)
    prfunc <- function(i) {
        res <- glmnet(features_train[tr_inds, ], as.numeric(train_resp[i, tr_inds]), standardize = FALSE)
        pr <- predict(res, features_train[te_inds, ], s=lambdas)
        pr
    }
    res <- mclapply(1:100, prfunc, mc.cores = 30)
                                        #res <- lapply(1:100, prfunc)
    pr_error <- numeric(nlambdas)
    misc_error <- numeric(nlambdas)
    for (i in 1:nlambdas) {
        pvalid <- matrix(0, nte, 100)
        for (j in 1:100) {
            pvalid[, j] <- res[[j]][, i]
        }
        pr_error[i] <- sum((t(pvalid) - train_resp[, te_inds])^2)
        te_cl <- knn(pvalid %*% omega_e, t(train_resp[, te_inds]) %*% omega_e, 1:nte, k=1)
        misc_error[i] <- misc_error[i] + sum(te_cl != 1:nte)
    }
    misc_errors[ii, ] <- misc_error
    pr_errors[ii, ] <- pr_error
    print(ii)
}
proc.time()

misc_error <- apply(misc_errors, 2, mean)
pr_error <- apply(pr_errors, 2, mean)

lambdas[order(misc_error)[1]]
lambdas[order(pr_error)[1]]

saveRDS(misc_errors, "misc_error.rds")
saveRDS(pr_errors, "pr_error.rds")

#plot(lambdas, misc_error)
#plot(lambdas, pr_error)
