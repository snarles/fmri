library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')

hp <- function(p) {
  res <- -p * log(p) - (1-p) * log(1-p)
  res[p==0 | p== 1] <- 0
  res
}

logist_ident <- function(Sigma, K, mc.reps = 1000) {
  p <- dim(Sigma)[1]
  mcs <- sapply(1:mc.reps,
                function(i) {
                  X <- mvrnorm(K, rep(0, p), Sigma)
                  ps <- 1/(1 + exp(-X))
                  Y <- (rand(K, p) < ps) + 0
                  nc <- rowSums(log(1 - ps))
                  pr <- t(t(Y %*% t(X)) + nc)
                  (lbls <- apply(pr, 1, function(v) order(-v)[1]))
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}


logist_mi <- function(Sigma, mc.reps = 1e5, h_est = h_jvhw) {
  p <- dim(Sigma)[1]
  X <- mvrnorm(mc.reps, rep(0, p), Sigma)
  ps <- 1/(1 + exp(-X))
  ces <- apply(ps, 1, function(v) {
    sum(hp(v))
  })
  Y <- (rand(mc.reps, p) < ps) + 0
  ylabs <- apply(Y, 1, function(v) {
    paste(v, collapse = "")
  })
  tab <- table(ylabs)
  #print(table(tab))
  #tab <- tab/sum(tab)
  #hy <- sum(-tab * log(tab))
  hy <- h_est(tab)
  hy - mean(ces)
}





