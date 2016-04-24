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

k.comps <- 10 ## number of components of mixture
ps <- rbeta(k.comps, 1, 1) ## Ps
prob.ps <-  rbeta(k.comps, 1, 1) ## prob of ps
prob.ps <- prob.ps/sum(prob.ps)

k <- 10 # size of binom
n <- 40 # how many dudes you see = k * r
Ys <- generate_sample(ps, prob.ps, k, n)
hist(ys)

pk_moment(ps, prob.ps, 2)
pk_moment(ps, prob.ps, 3)
