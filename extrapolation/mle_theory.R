####
##  Maximum likelihood deconvolution
####

generate_sample <- function(ps, prob.ps, k, n) {
  samp.p <- sample(ps, n, TRUE, prob = prob.ps)
  obs <- rbinom(n, k, samp.p)
  obs
}

pk_moment <- function(ps, prob.ps, k) {
  sum(ps^k * prob.ps)
}

est_moment <- function(res, k) {
  sum(res$theta[, 2]^k * res$lambda)
}

k.comps <- 3 ## number of components of mixture
ps <- rbeta(k.comps, 1, 1) ## Ps
prob.ps <-  rbeta(k.comps, 1, 1) ## prob of ps
prob.ps <- prob.ps/sum(prob.ps)

k <- 20 # size of binom
n <- 1000 # how many dudes you see = k * r
Ys <- generate_sample(ps, prob.ps, k, n)
Ymat <- cbind(k - Ys, Ys)
hist(Ys/k)
library(mixtools)
res <- multmixEM(Ymat, k = 10, verb = TRUE)
est_moment(res, 2)
pk_moment(ps, prob.ps, 2)
mean(sapply(Ys, binmom, k, 2))

est_moment(res, 3)
pk_moment(ps, prob.ps, 3)
mean(sapply(Ys, binmom, k, 3))

est_moment(res, 5)
pk_moment(ps, prob.ps, 5)
mean(sapply(Ys, binmom, k, 5))

est_moment(res, 10)
pk_moment(ps, prob.ps, 10)
mean(sapply(Ys, binmom, k, 10))


est_moment(res, 20)
pk_moment(ps, prob.ps, 20)
mean(sapply(Ys, binmom, k, 20))

est_moment(res, 25)
pk_moment(ps, prob.ps, 25)

est_moment(res, 30)
pk_moment(ps, prob.ps, 30)
mean((Ys/k)^30)


mean((Ys/k)^2)

rbind(res$theta[, 2][order(-res$lambda)], sort(res$lambda, decreasing = TRUE))
rbind(ps, prob.ps)
