library(pracma)
library(MASS)
library(class)
library(parallel)
#library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')
#source("approximation/extreme_value.R")



binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

## log 1+x
log_1_plus <- function(x) {
  if (abs(x) > 1e-3) return(log(1 + x))
  sum(-1 * (-x)^(1:20)/(1:20))
}

get_sub_errs <- function(pmat, true_ys, ks) {
  k <- nrow(pmat)
  p2 <- apply(pmat, 2, rank)
  true_ranks <- p2[cbind(true_ys, 1:ncol(pmat))]
  1 - sapply(ks, function(v) mean(binmom(true_ranks - 1, k - 1, v - 1)))
}


## mu ~ N(0, I)
## y ~ N(mu*, sigma^2 I)

## most naive implementation of mc() for identity cov
mc_ident_fs <- function(p, sigma2, sigma2_tr, K, mc.reps = 1000) {
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- randn(K, p)
                  ys <- mus + sqrt(sigma2) * randn(K, p)
                  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
                  lbls <- knn(mu_hats, ys, cl=1:K)
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}

mc_ident_fs_curve <- function(p, sigma2, sigma2_tr, K, mc.reps = 100) {
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- randn(K, p)
                  ys <- mus + sqrt(sigma2) * randn(K, p)
                  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
                  pmat <- -pdist2(mu_hats, ys)
                  get_sub_errs(pmat, 1:K, 1:K)
                })
  rowMeans(mcs)
}

## uses noncentral chi squared
## same Evalue as mc_ident
mc_ident2 <- function(p, sigma2, K, mc.reps = 1000) {
  alpha <- sigma2/(1 + sigma2)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y2 <- (1 + sigma2) * rchisq(1, df = p)
                  d1 <- alpha * rchisq(1, df = p, ncp = alpha * y2)
                  ds <- rchisq(K - 1, df = p, ncp = y2)
                  (min(ds) < d1)
                })
  mean(mcs)
}

# ### Which is faster?
# 
# p <- 3
# sigma2 <- 0.1
# sigma2_tr <- 0.1
# K <- 500
# mc.reps <- 100
# 
# t1 <- proc.time()
# res1 <- numeric(K)
# for (k in 2:K) res1[k] <- mc_ident_fs(p, sigma2, sigma2_tr, k, mc.reps)
# proc.time() - t1
# 
# t1 <- proc.time()
# res2 <- mc_ident_fs_curve(p, sigma2, sigma2_tr, K, mc.reps)
# proc.time() - t1
# 
# plot(res1); points(res2, col= "red")
# 
# ### Which is faster? 2
# 
# t1 <- proc.time()
# res1 <- numeric(K)
# for (k in 2:K) res1[k] <- mc_ident2(p, sigma2, k, mc.reps)
# proc.time() - t1
# 
# 
# t1 <- proc.time()
# res2 <- mc_ident_fs_curve(p, sigma2, 0, K, mc.reps)
# proc.time() - t1
# 
# res1 <- numeric(K)
# for (k in 2:K) res1[k] <- mc_ident2(p, sigma2, k, mc.reps * 10)
# 
# plot(res1); points(res2, col= "red")

## fitting the piK function

mugrid <- seq(0.1, 5, 0.1)
sigma2grid <- seq(0.1, 5, 0.1)
length(mugrid) * length(sigma2grid)

parmat <- cbind(mu = rep(mugrid, each = length(sigma2grid)),
                sigma2 = rep(sigma2grid, length(mugrid)))

nrow(parmat)

t1 <- proc.time()
par_res <- mclapply(1:nrow(parmat), function(i) {
  mu <- parmat[i, 1]
  sigma2 <- parmat[i, 2]
  1-piK(K = 1:250, mc.reps = 1000, mu, sigma2)
}, mc.cores = 2)
proc.time() - t1

precomputed_curves <- do.call(cbind, par_res)
save(parmat, precomputed_curves, file = "approximation/ident_curve_precomp.rda")

K <- 250
k_sub <- 125
p <- 3
sigma2 <- 0.1
mc.reps <- 100
acs <- 1-mc_ident_fs_curve(p, sigma2, 0, K, mc.reps)

## find best 2-parameter fit
dd <- colSums((precomputed_curves - acs)[1:k_sub, ]^2)
(pars <- parmat[which.min(dd), ])
acs_hat_2p <- precomputed_curves[, which.min(dd)]

## find best fit with sigma=1
dd[parmat[, "sigma2"] != 1] <- Inf
mu_fit <- parmat[which.min(dd), 1]
acs_hat_1p <- precomputed_curves[, which.min(dd)]

plot(1:250, acs, type = "l", ylim = c(0, 1))
lines(acs_hat_2p, col = "blue")
lines(acs_hat_1p, col = "red")
abline(v = k_sub)
legend(150, 0.5, col = c("black", "red", "blue"), lwd = 1, 
       legend = c("empirical", "1par", "2par"))
