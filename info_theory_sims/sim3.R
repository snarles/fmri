####
##  Simulation with estimation
####

library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')
source("info_theory_sims/methods_source.R")
source("entropy/jvhw.R")
source("entropy/gastpar.R")
hp <- function(p) {
  res <- -p * log(p) - (1-p) * log(1-p)
  res[p==0 | p== 1] <- 0
  res
}

logit_class <- function(Y, ps) {
  logodds <- log(ps/(1-ps))
  nc <- rowSums(log(1 - ps))
  pr <- t(t(Y %*% t(logodds)) + nc)
  lbls <- apply(pr, 1, function(v) order(-v)[1])
  lbls
}

####
##  Simulation code
####

## arguments


run_simulation <- function(p, q, Bmat, 
                           data.reps, m.folds, k.each, r.each, r.train,
                           mc.reps = 1e5, h_est = h_mle) {
  r.test <- r.each - r.train
  ## compute true MI
  
  X <- randn(mc.reps, p)
  ps <- 1/(1 + exp(-X %*% Bmat))
  ces <- apply(ps, 1, function(v) {
    sum(hp(v))
  })
  Y <- (rand(mc.reps, q) < ps) + 0
  ylabs <- apply(Y, 1, function(v) {
    paste(v, collapse = "")
  })
  tab <- table(ylabs)
  hy <- h_est(tab)
  (mi_true <- hy - mean(ces))
  
  ## draw data
  
  res <- list()
  for (i in 1:data.reps) {
    Xall <- numeric(0)
    Yall <- numeric(0)
    mcs <- numeric(m.folds)
    mi_cms <- numeric(m.folds)
    
    for (j in 1:m.folds) {
      X <- randn(k.each, p)
      ps <- 1/(1 + exp(-X %*% Bmat))
      ## generate data
      zs <- rep(1:k.each, each = r.each)
      test_filt <- rep(c(rep(FALSE, r.train), rep(TRUE, r.test)), k.each)
      Y <- (rand(k.each * r.each, q) < ps[zs, ]) + 0
      mu.tr <- zeros(k.each, q)
      for (ii in 1:k.each) {
        mu.tr[ii, ] <- colMeans(Y[zs ==ii, ][1:r.train, ])
      }
      ## Bayes smoothing
      mu.tr <- mu.tr * (r.train)/(r.train + 1) + .5/(r.train+1)
      ## prediction
      Y.te <- Y[test_filt, ]
      zs.te <- zs[test_filt]
      lbls.bayes <- logit_class(Y.te, ps)
      lbls <- logit_class(Y.te, mu.tr)
      (mc <- sum(lbls != zs.te)/length(zs.te))
      cm <- table(zs.te, lbls); cm <- cm/sum(cm)
      (mi_cm <- cm_to_I(cm))
      ## accumulation
      mcs[j] <- mc
      mi_cms[j] <- mi_cm
      Xall <- rbind(Xall, X)
      Yall <- rbind(Yall, Y)
    }
    ## ML-based estimates
    (mi_cm <- mean(mi_cms))
    abe <- mean(mcs)
    (mi_fano <- Ihat_fano(abe, k.each))
    (mi_ls <- Ihat_LS(abe, k.each))
    ## nonparametric estimate of MI
    ylabs <- apply(Yall, 1, function(v) paste(v, collapse = ""))
    yids <- as.numeric(as.factor(ylabs))
    ytab <- t(matrix(yids, r.each, m.folds * k.each))
    (mi_0 <- mi_naive(ytab))
    (mi_5 <- anthropic_correction(ytab, 0.5))
    (mi_9 <- anthropic_correction(ytab, 0.9))
    (mi_j <- mi_naive(ytab, h_jvhw))
    ## accumulation
    res[[i]] <- c(mi_cm = mi_cm, mi_fano = mi_fano, mi_ls = mi_ls,
                  mi_0 = mi_0, mi_5 = mi_5, mi_9 = mi_9, mi_j = mi_j)
  }
  res <- do.call(rbind, res)
  cbind(mi_true, res)
}


####
##  DEMO
####

p <- 5; q <- 10
Bmat <- 3 * randn(p, q)
data.reps <- 3
m.folds <- 3
k.each <- 3
r.each <- 20
r.train <- floor(0.8 * r.each)
(N = m.folds * k.each * r.each)
mc.reps <- 1e5
#h_est <- h_jvhw
h_est <- h_mle
res <- run_simulation(p, q, Bmat, 
                      data.reps, m.folds, k.each, r.each, r.train, 
                      mc.reps, h_est)
res

