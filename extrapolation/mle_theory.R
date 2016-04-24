####
##  Maximum likelihood deconvolution
####

generate_sample <- function(ps, prob.ps, k, n) {
  samp.p <- sample(ps, n, TRUE, prob = ps)
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

k <- 100 # size of binom
n <- 1000 # how many dudes you see = k * r
Ys <- generate_sample(ps, prob.ps, k, n)
Ymat <- cbind(k - Ys, Ys)
hist(Ys)
library(mixtools)
res <- multmixEM(Ymat, rep(1/10, 10), theta = NULL)
est_moment(res, 2)
pk_moment(ps, prob.ps, 2)
est_moment(res, 3)
pk_moment(ps, prob.ps, 3)



