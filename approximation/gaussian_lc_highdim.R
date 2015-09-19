####
##  High-dimensional learning curves
####

source('approximation//gaussian_lc_source.R')
source('approximation//gaussian_2c_source.R')
mcc <- 7 ## number of mc.cores


# mean of exp(v)
meanexp <- function(v) {
  vm <- max(v)
  mean(exp(v - vm)) * exp(vm)
}

## High-dimensional asymptotic misclassification curve for Sigma = cc/sqrt(d) * eye(d)
mcK_I <- function(cc, maxK, mc.reps = 1e4) {
  samp <- qnorm(((1:mc.reps) - 0.5)/mc.reps)
  temp <- log(1 - pnorm(samp - sqrt(cc)))
  tmat <- ((1:maxK) - 1) * repmat(temp, maxK, 1)
  1 - apply(tmat, 1, meanexp)
}

sdS <- function(cc, d) sqrt(2) * cc/sqrt(d)
#SigHd <- function(cc, d) cc/d * eye(d)
SigHd <- function(cc, d) {
  Sigma <- cov(randn(2*d, d))
  cc/sum(diag(Sigma)) * Sigma
}


d <- 100
cc <- 0.1
sdS(cc, d)

mcK_I(cc, 5)
mcK(cc/d * eye(d), 5)
mc2c_3(cc/d * eye(d))
Sigma <- SigHd(cc, d)

maxK <- 20
scalings <- floor((1:20)^1.5)
Sigmas <- lapply(scalings, function(s) SigHd(cc, s))
res <- mclapply(Sigmas, function(s) mcK(s, maxK), mc.cores = mcc)
mat <- do.call(cbind, res)
par(bg = "grey")
matplot(mat, type = "l", col = rainbow(length(Sigmas)), lty = 1, lwd = 1)
lines(1:maxK, mcK_I(cc, maxK), lwd = 2)

