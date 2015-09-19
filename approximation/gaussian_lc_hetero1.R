####
##  DIVERSITY of LARGE K/LARGE TR/SMALL D case
##  Do curves vary along only two parameters? (tr and dim)?
####

source('approximation//gaussian_lc_source.R')
source('approximation//gaussian_2c_source.R')
mcc <- 7
par(bg = "grey")

eps <- exp((-10:0)/3)
Sigma0 <- diag(c(1, 1, 0, 0))
Sigma1s <- lapply(eps, function(ep) diag(c(1, ep, ep, ep)))

c0 <- build_mc2c_table(Sigma0, goal=1/4)
Sigma0 <- c0 * Sigma0
mc2c_3(Sigma0)
mcK(Sigma0, 10)

c1s <- sapply(Sigma1s, function(s) build_mc2c_table(s, goal=1/4))
Sigma1s <- lapply(1:length(Sigma1s), function(i) c1s[i] * Sigma1s[[i]])
sapply(Sigma1s, mc2c_3)
mcK(Sigma1s[[1]], 10)

maxK <- 10
res <- mclapply(Sigma1s, function(s) mcK(s, maxK, N = 200, mc.reps=1000), mc.cores = mcc)
mat <- do.call(cbind, res)
matplot(mat, type = "l", col = rainbow(length(Sigmas)), lty = 1, lwd = 1)
lines(1:maxK, mcK(Sigma0, maxK, N = 200, mc.reps = 1000), lwd = 2)

## need better approximation of mc2c!
