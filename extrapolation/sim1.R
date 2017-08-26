source("extrapolation/gaussian_source.R")
source("extrapolation/psuedolikelihood.R")

k <- 100
sigma <- 0.4
p <- 3
g.res <- 200
mus <- randn(k, p)
Ys <- get_vs(mus, sigma)
mean(binmom(Ys, k, k-1))
res_gu <- compute_gu(p, sigma, g.res = 1000)
gu0 <- res_gu$gu/sum(res_gu$gu)
us <- res_gu$mids
plot(us, gu0, type = "l")
sum(res_gu$mids ^ (k - 1) * res_gu$gu)/sum(res_gu$gu)
est <- fit_pm_models(Ys, k, us = res_gu$mids)
lines(res_gu$mids, est$gu_mple, col = "red")
lines(res_gu$mids, est$gu_mom, col = "blue")
lines(res_gu$mids, est$gu_mono, col = "green")
lines(res_gu$mids, est$gu_mm, col = "orange")

K <- 500
mus2 <- rbind(mus, randn(K - k, p))
(accK <- get_bayes_acc(mus2, sigma))
estimates <- c(accK = accK,
               accBayes = sum(us^(K-1) * gu0),
               accMPLE = sum(us^(K-1) * est$gu_mple),
               accMOM = sum(us^(K-1) * est$gu_mom),
               accMM = sum(us^(K-1) * est$gu_mm))
estimates

