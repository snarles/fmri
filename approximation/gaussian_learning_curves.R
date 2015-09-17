####
##  Learning curves for multi-class gaussian classification
####

source('approximation//gaussian_lc_source.R')
source('approximation//gaussian_2c_source.R')
mcc <- 7 ## number of mc.cores

####
##  Some initial experiments
##  to see how well MI correlates with classification (not well)
####

Kmax <- 20
N <- 200
mc.reps <- 1000
.5/sqrt(mc.reps)

expment1 <- function(seed) {
  set.seed(seed)
  Sigma <- rchisq(1, 5) * cov(randn(8, 5))
  curve <- mcK(Sigma, Kmax, N, mc.reps)
  mutual <- mi(Sigma)
  list(curve = curve, mutual = mutual, Sigma = Sigma)
}

expment2 <- function(seed) {
  set.seed(seed)
  Sigma <- scalings[seed] * Sigma0
  curve <- mcK(Sigma, Kmax, N, mc.reps)
  mutual <- mi(Sigma)
  list(curve = curve, mutual = mutual, Sigma = Sigma)
}

Sigma0 <- cov(randn(8, 5))
scalings <- (5:15)/10


res <- lclapply(1:100, expment1, mc.cores = 7)
#res <- lclapply(1:length(scalings), expment2, mc.cores = 7)
zattach(res)
curves <- do.call(cbind, curve)
mutual <- unlist(mutual)
cols <- rainbow(length(mutual))[order(mutual)]
par(bg = "grey")
plot(1:Kmax, 1:Kmax/Kmax, pch = ".", col = "grey")
for (ii in 1:length(Sigma)) {
  lines(1:Kmax, curve[[ii]], col = cols[[ii]])
}


trs <- sapply(Sigma, function(s) sum(diag(s)))
tr2s <- sapply(Sigma, function(s) sum(diag(s^2)))
res <- lm(log(curves[2, ]) ~ trs + tr2s + trs^2 + tr2s^2)
plot(log(curves[2, ]), res$fitted.values)

plot(trs, log(curves[2, ]))
plot(tr2s, log(curves[2, ]))

plot(mutual, log(curves[2, ]))
plot(curves[2, ], curves[10, ])
plot(curves[2, ], curves[20, ])
plot(curves[10, ], curves[20, ])


par(bg = "white")


####
##  Experiments to see how well
##  two-class performance is correlated to K-class performance
####

####
##  Part 1: Effect of dimensionality
####

## scalings chosen to have mc(., 2) = 1/4
scalings <- unlist(mclapply(1:20,
                            function(k) build_mc2c_table(eye(k), 1/4), mc.cores = mcc))
plot(scalings * 1:20)
res <- mclapply(1:20, function(s) mcK(scalings[s] * eye(s), 40), mc.cores = 3)
mat <- do.call(cbind, res)
par(bg = "grey")
matplot(1:40, mat, type = "l", col = rainbow(20), lty = 1, lwd = 1)
lines(1:40, 1 - 1/(1:40), lwd = 2)

plot(mat[20, ])

#### 
##  Part II: Interpolating between dimensionality
####


fracs <- c(0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, .5)
fracs <- c(fracs, 1 - rev(fracs)[-1])
#fracs <- 0:20/20
Sigma0s <- lapply(fracs, function(x) x * matA + (1-x) * matB)
scalings <- unlist(mclapply(Sigma0s,
                            function(s) build_mc2c_table(s, 1/4),
                            mc.cores = mcc))
Sigmas <- lapply(1:length(Sigma0s), function(i) scalings[i] * Sigma0s[[i]])
res <- mclapply(Sigmas, function(s) mcK(s, 40), mc.cores = mcc)
mat <- do.call(cbind, res)
par(bg = "grey")
matplot(mat, type = "l", col = rainbow(length(fracs)), lty = 1, lwd = 1)
plot(fracs, mat[40, ])
fracs
lapply(Sigmas, diag)
