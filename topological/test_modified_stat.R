####
##  Try unbiased stat S
####

source("topological/rsa_boot_source.R")

p <- 5
q <- 2
SigmaX <- 3 * cov(randn(2*q, q)); SigmaY <- 3 * cov(randn(2 * q, q))
SigmaE1 <- 10 * cov(randn(5 * p, p)); SigmaE2 <- cov(randn(5 * p, p))
A_0 <- randn(p, q); B_0 <- svd(randn(p, p))$u %*% A_0
h0_small <- regression_data_model_(A_0, B_0, SigmaX, SigmaY, SigmaE1, SigmaE2)

B_1 <- randn(p, q)
h1_small <- regression_data_model_(A_0, B_1, SigmaX, SigmaY, SigmaE1, SigmaE2)

dat <- h0_small(200, 200)
mus <- sample_moments(dat)
c(f2(mus$Ahat, A_0), f2(mus$Bhat, B_0))


nX <- 100; nY <- 100; mc.reps = 1000
res0 <- h0_small(nX, nY)
res1 <- h1_small(nX, nY)

c(inverse_bca_test(res0, stat.T, mc.reps), inverse_bca_test(res1, stat.T, mc.reps))
c(inverse_bca_test(res0, stat.S, mc.reps), inverse_bca_test(res1, stat.S, mc.reps))
c(inverse_bca_test(res0, stat.Su, mc.reps), inverse_bca_test(res1, stat.Su, mc.reps))
