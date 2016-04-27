####
## factor model for generating probs
####

source("extrapolation/example2d.R")
source("extrapolation/mle_theory.R")
source("extrapolation/mcmc.R")
source("extrapolation/constrained_mle.R")

library(pracma)

nX <- 500
nY <- 500
rk <- 10
xfax <- randn(nX, rk)
yfax <- randn(nY, rk)
ipmat <- xfax %*% t(yfax)
coefs <- c(-2, 1, 0, 0)
pmat <- exp(coefs[1] + coefs[2] * ipmat + coefs[3] * ipmat^2 + coefs[4] * ipmat)
for (i in 1:10) {
  rs <- rowSums(pmat)
  pmat <- pmat/rs
  pmat <- t(pmat)
}
# avg_mc_acc_naive(pmat, 20)
# avg_mc_acc_p(pmat, 20)
# avg_mc_acc_naive(pmat, 40)
# avg_mc_acc_p(pmat, 40)
# avg_mc_acc_p(pmat, 90)
# avg_mc_acc_p(pmat, 300)


####
##  LOOKING AT TRUE DISTIRUBITION
####

rankconv <- (apply(pmat, 2, rank) - 0.5)/nrow(pmat)
plot(density(as.numeric(rankconv), weights = as.numeric(pmat)/sum(pmat)))
ps <- sample(as.numeric(rankconv), 10000, replace = TRUE, prob = as.numeric(pmat)/sum(pmat))
plot(sort(ps), type = "l")
saveRDS(ps, file = "extrapolation/testcase_ps.rds")

####
##  USING MLE
####

ppmat <- empirical_p_dist(pmat, 30, 100)
table(as.numeric(ppmat))
res <- res_mixtools(ppmat, 30)
pseq <- seq(0, 1, 1/300)
cm <- cons_mle_est(ppmat, k = 30, ps = pseq, lbda = 0.1)

avg_mc_acc_p(pmat, 30)
est_moment(res, 30)
mean(sapply(as.numeric(ppmat), binmom, k, 30))
cm_est_moment(cm, 30)

hist(ps)
plot(cm$gu)

avg_mc_acc_p(pmat, 300)
est_moment(res, 300)

####
##  USING MCMC
# ####
# 
# mc <- mcmc_fitting(ppmat, 30, iter = 50, chains = 2)
# mean(mchain_est(mc, 90))
# avg_mc_acc_p(pmat, 90)
# 
# mean(mchain_est(res, 300))
# avg_mc_acc_p(pmat, 300)

