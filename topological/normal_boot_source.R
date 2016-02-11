####
##  Bootstrapping testing of procrustes problem
##  Generates A, B matrices with general covariance patterns
####

source("rsa_boot_source.R")
## DATA MODEL

n_data_model_ <- function(A, B, SigmaX, SigmaY) {
  p <- dim(A)[1]; q <- dim(A)[2]
  sampler <- function(nX, nY) {
    noise_mult <- 1
    if (nX == 0 && nY == 0) {
      nX <- 1; nY <- 1
      noise_mult <- 0
    }
    rawX <- mvrnorm(nX, as.numeric(A), SigmaX)
    rawY <- mvrnorm(nY, as.numeric(B), SigmaY)
    dat <- rbind(cbind(0, rawX), cbind(1, rawY))
    list(p = p, q  = q, dat = dat)
  }
  # force eval
  lala <- sampler(10, 10)
  sampler
}

n_sample_moments <- function(res) {
  p <- res$p; q <- res$q; dat <- res$dat
  rawX <- dat[dat[,1] == 0, -1, drop = FALSE]
  rawY <- dat[dat[,1] == 1, -1, drop = FALSE]
  Ahat <- matrix(colMeans(rawX), p, q)
  Bhat <- matrix(colMeans(rawY), p, q)
  list(Ahat = Ahat, Bhat = Bhat)
}

## TEST STATISTICS

n_stat.T <- function(res) {
  mus <- n_sample_moments(res)
  res <- svd(mus$Ahat %*% t(mus$Bhat))
  R <- res$v %*% t(res$u)
  stat.T.raw <- R %*% mus$Ahat - mus$Bhat
  as.numeric(stat.T.raw)
}

n_stat.S <- function(res) {
  mus <- n_sample_moments(res)
  stat.S.raw <- t(mus$Ahat) %*% mus$Ahat - t(mus$Bhat) %*% mus$Bhat
  as.numeric(stat.S.raw)
}

####
##  Demo
####

## REGRESSION

p <- 2
q <- 5
SigmaX <- 3 * cov(randn(2*q, q)); SigmaY <- 3 * cov(randn(2 * q, q))
A_0 <- randn(p, q); B_0 <- svd(randn(p, p))$u %*% A_0
h0_small <- regression_data_model_(A_0, B_0, SigmaX, SigmaY)

B_1 <- randn(p, q)
h1_small <- regression_data_model_(A_0, B_1, SigmaX, SigmaY)

dat <- h0_small(200, 200)
mus <- sample_moments(dat)
c(f2(mus$Ahat, A_0), f2(mus$Bhat, B_0))


nX <- 100; nY <- 100; mc.reps = 1000
res0 <- h0_small(nX, nY)
res1 <- h1_small(nX, nY)

c(inverse_bca_test(res0, stat.T, mc.reps), inverse_bca_test(res1, stat.T, mc.reps))
c(inverse_bca_test(res0, stat.S, mc.reps), inverse_bca_test(res1, stat.S, mc.reps))

## NORMAL

p <- 2
q <- 30
SigmaX <- 3 * cov(randn(2*q, q)); SigmaY <- 3 * cov(randn(2 * q, q))
SigmaE <- cov(randn(2 * p, p))
SigmaA <- SigmaE %x% SigmaX
SigmaB <- SigmaE %x% SigmaY
A_0 <- randn(p, q); B_0 <- svd(randn(p, p))$u %*% A_0
h0_small <- n_data_model_(A_0, B_0, SigmaA, SigmaB)

B_1 <- randn(p, q)
h1_small <- n_data_model_(A_0, B_1, SigmaA, SigmaB)

dat <- h0_small(200, 200)
mus <- n_sample_moments(dat)
c(f2(mus$Ahat, A_0), f2(mus$Bhat, B_0))


nX <- 100; nY <- 100; mc.reps = 1000
res0 <- h0_small(nX, nY)
res1 <- h1_small(nX, nY)

c(inverse_bca_test(res0, n_stat.T, mc.reps), inverse_bca_test(res1, n_stat.T, mc.reps))
c(inverse_bca_test(res0, n_stat.S, mc.reps), inverse_bca_test(res1, n_stat.S, mc.reps))
