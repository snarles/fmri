####
##  High-dimensional learning curves
####

source('approximation//gaussian_lc_source.R')
source('approximation//gaussian_2c_source.R')
mcc <- 7 ## number of mc.cores




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


####
##  Increase C and number of classes
####


## C = 2 is fine w/ > 100 dimensions
cc <- 2
plot(mcK_I(cc, 50), type = "l")
maxK <- 50
scalings <- floor((1:20)^1.5)
Sigmas <- lapply(scalings, function(s) SigHd(cc, s))
res <- mclapply(Sigmas, function(s) mcK(s, maxK, N = 2000, mc.reps=10), mc.cores = mcc)
mat <- do.call(cbind, res)
par(bg = "grey")
matplot(mat, type = "l", col = rainbow(length(Sigmas)), lty = 1, lwd = 1)
lines(1:maxK, mcK_I(cc, maxK), lwd = 2)

## C = 3 shows a small gap
cc <- 3.5
plot(mcK_I(cc, 70), type = "l")
maxK <- 70
scalings <- floor((1:20)^1.5)
Sigmas <- lapply(scalings, function(s) SigHd(cc, s))
res <- mclapply(Sigmas, function(s) mcK(s, maxK, N = 2000, mc.reps=10), mc.cores = mcc)
mat <- do.call(cbind, res)
par(bg = "grey")
matplot(mat, type = "l", col = rainbow(length(Sigmas)), lty = 1, lwd = 1)
lines(1:maxK, mcK_I(cc, maxK), lwd = 2)


## C = 5 shows a gap at 100 dims
cc <- 5
plot(mcK_I(cc, 200), type = "l")
maxK <- 200
scalings <- floor((1:20)^1.5)
Sigmas <- lapply(scalings, function(s) SigHd(cc, s))
res <- mclapply(Sigmas, function(s) mcK(s, maxK, N = 2000, mc.reps=10), mc.cores = mcc)
mat <- do.call(cbind, res)
par(bg = "grey")
matplot(mat, type = "l", col = rainbow(length(Sigmas)), lty = 1, lwd = 1)
lines(1:maxK, mcK_I(cc, maxK), lwd = 2)

## C = 10 already a massive gap
cc <- 10
plot(mcK_I(cc, 1000), type = "l")
maxK <- 1000
scalings <- floor((1:20)^1.5)
Sigmas <- lapply(scalings, function(s) SigHd(cc, s))
res <- mclapply(Sigmas, function(s) mcK(s, maxK, N = 2000, mc.reps=10), mc.cores = mcc)
mat <- do.call(cbind, res)
par(bg = "grey")
matplot(mat, type = "l", col = rainbow(length(Sigmas)), lty = 1, lwd = 1)
lines(1:maxK, mcK_I(cc, maxK), lwd = 2)

matplot(mat[1:10, ], type = "l", col = rainbow(length(Sigmas)), lty = 1, lwd = 1)
lines(1:10, mcK_I(cc, maxK)[1:10], lwd = 2)

