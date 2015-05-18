#############################
## REDUCED RANK REGRESSION ##
#############################

## PREPARATION ##

library(magrittr)
library(Rcpp)
library(parallel)
library(glmnet)

isqrtm <- function(m) {
  res <- eigen(m)
  d <- res$values
  if (min(d) < -1e-5) warning("Negative eigenvalues in isqrtm")
  d[d < 0] <- 0
  d[d > 0] <- 1/sqrt(d[d > 0])
  v <- res$vectors
  return (v %*% diag(d) %*% t(v))
}


sourceCpp('ident_regression/data/pdist.cpp')

ddir <- "~/stat312data"
list.files(ddir)

load(paste0(ddir, "/feature_train.RData"))
load(paste0(ddir, "/train_resp_all.RData"))

complete_filt <- train_resp %>% colSums %>% is.na %>% `!`

resp <- train_resp[, complete_filt]
dim(resp) # 1750 22733

n_vx <- dim(resp)[2]

#res <- svd(feature_train)
#res$d

## SET UP SPLIT and also subsample voxels ##

n_svx <- 3000
set.seed(0)
train_inds <- sample(1750, 1600)
vox_inds <- sample(n_vx, n_svx)
tr_y <- resp[train_inds, vox_inds]
te_y <- resp[-train_inds, vox_inds]
tr_x <- feature_train[train_inds, ]
te_x <- feature_train[-train_inds, ]

svdx <- svd(cbind(1, tr_x))
regmat <- svdx$v %*% diag(svdx$d^(-2)) %*% t(svdx$v) %*% t(cbind(1, tr_x)) 
dim(regmat)
## REGULAR REGRESSION ##

i <- 1
el_err <- function(i) {
  res <- cv.glmnet(tr_x, tr_y[, i], alpha = .5)
  pr_y <- predict(res, te_x)
  err <- sum((te_y[, i] - pr_y)^2)
  err  
}

proc.time()
errs <- mclapply(1:50, el_err)
proc.time()
unlist(errs)

## REDUCED RANK REG ##

res_s <- svd(tr_y)
res_s$d

rho <- 0.2 # exponent for adaptive NN
lbda <- 1500 # lambda for adaptive NN
d2 <- pmax(0, res_s$d - lbda * res_s$d ^ (-rho))
tr_y_rr <- res_s$u %*% diag(d2) %*% t(res_s$v)


el_err2 <- function(i) {
  cff <- regmat %*% tr_y_rr[, i]
  pr_y <- cbind(1, te_x) %*% cff
  err <- sum((te_y[, i] - pr_y)^2)
  err  
}
el_err2(2)

proc.time()
errs2 <- mclapply(1:50, el_err2)
proc.time()


## compare

plot(unlist(errs), unlist(errs2))
abline(0, 1)
