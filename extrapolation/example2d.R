###
## x, y uniform marginals
## classification error
## P-function
## extrapolation
####

#'
#' nres <- 1000  ## much larger than K--the target extrapolation
#' apar <- 1; bpar <- 2 ## parameters for Beta distribution
create_pmat <- function(nres, apar, bar, dstrength, nits = 10) {
  pmat <- matrix(rbeta(nres^2, apar, bpar), nrow = nres)
  pmat <- pmat + nres * dstrength * pracma::eye(nres)
  ## row and column standardize
  for (i in 1:nits) {
    rs <- rowSums(pmat)
    pmat <- pmat/rs
    pmat <- t(pmat)
  }
  pmat
}

## average mc acc
avg_mc_acc_naive <- function(pmat, k, reso = 1e3) {
  acs <- numeric(reso)
  for (i in 1:reso){
    xs <- sample(nrow(pmat), k, TRUE)
    ys <- numeric(k)
    for (j in 1:k) ys[j] <- sample(1:ncol(pmat), 1, prob = pmat[xs[j], ])
    subps <- pmat[xs, ys]
    idx <- apply(subps, 2, function(v) which(v == max(v))[1])
    acs[i] <- sum(idx == 1:k)/k
  }
  mean(acs)
}

## avg mc error using p-distribution
avg_mc_acc_p <- function(pmat, k) {
  rankconv <- (apply(pmat, 2, rank) - 0.5)/nrow(pmat)
  sum(pmat/sum(pmat) * rankconv^(k - 1))
}

## get an empirical P--distribution
## r = repeats per k
empirical_p_dist <- function(pmat, k, r) {
  xs <- sample(nrow(pmat), k, TRUE)
  ys <- matrix(0, k, r)
  subps <- pmat[xs, ]
  Pmat <- matrix(0, k, r)
  for (j in 1:k) {
    ys[j, ] <- sample(1:ncol(pmat), r, prob = pmat[xs[j], ], replace = TRUE)
    tmat <- apply(subps[, ys[j, ]], 2, rank, ties.method = "random")
    Pmat[j, ] <- tmat[j, ]
  }
  Pmat
}

## binomial moment fomula
binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

estimate_mc_acc <- function(pmat, ksamp, r, k, nreps = 1) {
  ests <- numeric(nreps)
  for (i in 1:nreps) {
    PPmat <- empirical_p_dist(pmat, ksamp, r)
    ests[i] <- mean(binmom(PPmat, ksamp, k))
  }
  ests
}

####
##  TESTS
####
# pmat <- create_pmat(2000, 2, 10, 0)
# rankconv <- (apply(pmat, 2, rank) - 0.5)/nrow(pmat)
# plot(density(as.numeric(rankconv), weights = as.numeric(pmat)/sum(pmat)))
# ps <- sample(as.numeric(rankconv), 10000, replace = TRUE, prob = as.numeric(pmat)/sum(pmat))
# plot(sort(ps), type = "l")
# 
# avg_mc_acc_naive(pmat, 2)
# avg_mc_acc_p(pmat, 2)
# avg_mc_acc_naive(pmat, 3)
# avg_mc_acc_p(pmat, 3)
# 
# k <- 4
# avg_mc_acc_p(pmat, k)
# ests <- estimate_mc_acc(pmat, 10, 8, k, 100)
# mean(ests)
# 


