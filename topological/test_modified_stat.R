####
##  Try unbiased stat S
####

source("topological/rsa_boot_source.R")

p <- 5
q <- 2
SigmaX <- 3 * cov(randn(2*q, q)); SigmaY <- 3 * cov(randn(2 * q, q))
SigmaE1 <- 10 * cov(randn(5 * p, p)); SigmaE2 <- cov(randn(5 * p, p))
A_0 <- randn(p, q); B_0 <- svd(randn(p, p))$u %*% A_0
MA_0 <- t(A_0) %*% A_0
MB_0 <- t(B_0) %*% B_0
MA_0 - MB_0
h0_small <- regression_data_model_(A_0, B_0, SigmaX, SigmaY, SigmaE1, SigmaE2)

stat.MA <- function(res) as.numeric(sample_moments(res)$M_A)

###
#  Is M_A hat unbiased?
###

ss <- sampling_dist(h0_small, stat.MA, 200, 200, mc.reps = 1000, samples = TRUE)
rowMeans(ss)
as.numeric(MA_0)
layout(matrix(1:4, 2, 2))
for (i in 1:4) {
  hist(ss[i, ] - MA_0[i])
}

ss <- sampling_dist(h0_small, stat.Su, 200, 200, mc.reps = 1000, samples = TRUE)
rowMeans(ss)
layout(matrix(1:4, 2, 2))
for (i in 1:4) {
  hist(ss[i, ] - MA_0[i])
}

###
#  Testing p-values
###

B_1 <- randn(p, q)
h1_small <- regression_data_model_(A_0, B_1, SigmaX, SigmaY, SigmaE1, SigmaE2)

nX <- 200; nY <- 200; mc.reps = 1000
res0 <- h0_small(nX, nY)
res1 <- h1_small(nX, nY)

c(inverse_bca_test(res0, stat.T, mc.reps), inverse_bca_test(res1, stat.T, mc.reps))
c(inverse_bca_test(res0, stat.S, mc.reps), inverse_bca_test(res1, stat.S, mc.reps))
c(inverse_bca_test(res0, stat.Su, mc.reps), inverse_bca_test(res1, stat.Su, mc.reps))

null_res <- bootstrap_sampling_dist(res0, stat.S, mc.reps)
layout(matrix(1:4, 2, 2)); for (i in 1:4) hist(null_res[i, ])

null_res <- bootstrap_sampling_dist(res0, stat.Su, mc.reps)
layout(matrix(1:4, 2, 2)); for (i in 1:4) hist(null_res[i, ])


