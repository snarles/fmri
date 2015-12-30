####
##  Simulation with estimation
####


library(pracma)
library(MASS)
library(class)
library(parallel)
#library(reginference) ## see github.com/snarles/misc
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

mi_ident_case <- function(p, mult, mc.reps = 1e3) {
  xs <- 3 * ((-mc.reps):mc.reps)/mc.reps * sqrt(2 * log(mc.reps))
  delt <- xs[2] - xs[1]
  ce <- delt * sum(dnorm(xs) * (log(1 + exp(-xs * mult))/(1 + exp(-xs * mult)) + 
                                  (log(1 + exp(xs * mult))/(1 + exp(xs * mult)))))
  p * (log(2) - ce)
}

####
##  Simulation code
####

## arguments


run_simulation <- function(Bmat, m.folds, k.each, r.each, r.train) {
  p <- dim(Bmat)[1]; q <- dim(Bmat)[2]
  r.test <- r.each - r.train
  ## draw data
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
    for (ii in 1:k.each) mu.tr[ii, ] <- colMeans(Y[zs ==ii, ][1:r.train, ])
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
  abe0 <- mean(mcs)
  n_te <- r.test * k.each * m.folds
  abe <- n_te/(n_te+1) * abe0 + (1 - 1/k.each)/(n_te + 1) # bayes smoothing
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
  ## answer
  c(mi_cm = mi_cm, mi_fano = mi_fano, mi_ls = mi_ls,
    mi_0 = mi_0, mi_5 = mi_5, mi_9 = mi_9, mi_j = mi_j, abe0 = abe0, abe = abe)
}

get_abe <- function(Bmat, k.each, r.each = 1, mc.reps=1e3, mcc=0) {
  p <- dim(Bmat)[1]; q <- dim(Bmat)[2]
  mcs <- mclapply0(1:mc.reps, function(i) {
    X <- randn(k.each, p)
    ps <- 1/(1 + exp(-X %*% Bmat))
    ## generate data
    zs <- rep(1:k.each, each = r.each)
    Y <- (rand(k.each * r.each, q) < ps[zs, ]) + 0
    lbls <- logit_class(Y, ps)
    sum(lbls != zs)/length(zs)
  }, mc.cores = mcc)
  abe <- mean(unlist(mcs))
  sd_abe <- sd(unlist(mcs))
  mc_b_ls <- Ihat_LS(abe, k.each)
  c(abe = abe, mc_b_ls = mc_b_ls, sd_abe = sd_abe)
}

## parallelize computing mi efficiently
compute_mi <- function(Bmat, mc.reps = 1e5, mcc = 0) {
  
  ces <- unlist(mclapply0(1:mcc, function(i) {
    X <- randn(mc.reps, p)
    ps <- 1/(1 + exp(-X %*% Bmat))
    apply(ps, 1, function(v) sum(hp(v)))
  }, mc.cores = mcc))
  
  ylabs <- unlist(mclapply0(1:mcc, function(i) {
    X <- randn(mc.reps, p)
    ps <- 1/(1 + exp(-X %*% Bmat))
    Y <- (rand(mc.reps, q) < ps) + 0
    apply(Y, 1, function(v) paste(v, collapse = ""))
  }, mc.cores = mcc))
  
  tab <- table(ylabs)
  c(mi_mle = h_mle(tab) - mean(ces), mi_jvhw = h_jvhw(tab) - mean(ces), ces = mean(ces))  
}

run_simulations <- function(Bmat, m.folds, k.each, r.each, r.train,
                            mcc = 0, data.reps = 10) {
  res <- mclapply0(1:data.reps, 
                   function(i) run_simulation(Bmat, m.folds, k.each, r.each, r.train),
                   mc.cores = mcc)
  do.call(rbind, res)
}





run_instance <- function(X, Yall, zs, nsub, ntr) {
  k.each <- dim(X)[1]
  Y.tr <- Yall[1:ntr, ]
  Y.te <- Yall[(ntr + 1):nsub, ]
  Y <- Yall[1:nsub, ]
  zs.sub <- zs[1:nsub]
  zs.tr <- zs[1:ntr]
  zs.te <- zs[(ntr + 1):nsub]
  mu.tr <- zeros(k.each, q)
  r.train <- ntr/k.each
  for (ii in 1:k.each) mu.tr[ii, ] <- colMeans(Y.tr[zs.tr ==ii, , drop = FALSE])
  ## Bayes smoothing
  mu.tr <- mu.tr * (r.train)/(r.train + 1) + .5/(r.train+1)
  ## prediction
  lbls <- logit_class(Y.te, mu.tr)
  (mc <- sum(lbls != zs.te)/length(zs.te))
  cm <- table(zs.te, lbls); cm <- cm/sum(cm)
  (mi_cm <- cm_to_I(cm))
  abe0 <- mc
  n_te <- nsub - ntr
  abe <- n_te/(n_te+1) * abe0 + (1 - 1/k.each)/(n_te + 1) # bayes smoothing
  (mi_fano <- Ihat_fano(abe, k.each))
  (mi_ls <- Ihat_LS(abe, k.each))
  ## nonparametric estimate of MI
  ylabs <- apply(Y, 1, function(v) paste(v, collapse = ""))
  yids <- as.numeric(as.factor(ylabs))
  ytab <- list()
  for (i in 1:k.each) ytab[[i]] <- yids[zs.sub == i]
  ytab <- do.call(rbind, ytab)
  (mi_0 <- mi_naive(ytab))
  (mi_5 <- anthropic_correction(ytab, 0.5))
  (mi_9 <- anthropic_correction(ytab, 0.9))
  (mi_j <- mi_naive(ytab, h_jvhw))
  ## answer
  c(mi_cm = mi_cm, mi_fano = mi_fano, mi_ls = mi_ls,
    mi_0 = mi_0, mi_5 = mi_5, mi_9 = mi_9, mi_j = mi_j, abe0 = abe0, abe = abe)
}


####
##  Logistic entropy calculation
####

eta <- function(x) log(1 + exp(x))
deta <- function(x) 1/(1 + exp(-x))
d2eta <- function(x) exp(-x)/(1 + exp(-x))^2

nll <- function(y, x, Bmat) {
  sum(x^2)/2 + sum(eta(-y * (t(Bmat) %*% x)))
}

opt_nll <- function(y, Bmat) {
  p <- dim(Bmat)[1]
  res <- optim(rep(0, p), function(x) nll(y, x, Bmat))
  res$par
}

pr_laplace <- function(y, Bmat, log.p = FALSE) {
  p <- dim(Bmat)[1]; q <- dim(Bmat)[2]
  Bt <- Bmat %*% diag(-y)
  x0 <- opt_nll(y, Bmat)
  mu <- as.numeric(t(Bt) %*% x0)
  dn <- deta(mu); d2n <- d2eta(mu)
  (l0 <- nll(y, x0, Bmat))
  #   gd <- x0 + Bt %*% dn
  hs <- eye(p) + Bmat %*% diag(d2n) %*% t(Bmat)
  #   as.numeric((1/sqrt(2*pi))^p * 
  #                exp(-l0 + 1/2 * t(gd) %*% solve(hs, gd)) * 
  #                sqrt(det(2 * pi * solve(hs))))
  if (log.p) return(-l0 - log(det(hs))/2)
  exp(-l0)/sqrt(det(hs))
}

compute_mi2 <- function(Bmat, mc.reps = 1e5, mcc = 0) {
  
  diffs <- unlist(mclapply0(1:mcc, function(i) {
    X <- randn(mc.reps, p)
    ps <- 1/(1 + exp(-X %*% Bmat))
    ces <- apply(ps, 1, function(v) sum(hp(v)))
    Y <- (rand(mc.reps, q) < ps) + 0
    logps <- apply(2 * Y - 1, 1, function(v) pr_laplace(v, Bmat, TRUE))
    -ces - logps
  }, mc.cores = mcc))
  
  c(mu = mean(diffs), se = sd(diffs)/sqrt(length(diffs)))
}

## pr laplace with importance sampling correction
pr_laplace2 <- function(y, Bmat, log.p = FALSE, in.reps = 1e2) {
  p <- dim(Bmat)[1]; q <- dim(Bmat)[2]
  Bt <- Bmat %*% diag(-y)
  x0 <- opt_nll(y, Bmat)
  mu <- as.numeric(t(Bt) %*% x0)
  dn <- deta(mu); d2n <- d2eta(mu)
  (l0 <- nll(y, x0, Bmat))
  hs <- eye(p) + Bmat %*% diag(d2n) %*% t(Bmat)
  ## compute correction factor
  xsamp <- mvrnorm(in.reps, x0, solve(hs))
  #nlls <- apply(xsamp, 1, function(x) nll(y, x, Bmat))
  bx <- xsamp %*% Bt
  nlls <- rowSums(xsamp^2)/2 + rowSums(eta(bx))
  #prox_nlls <- apply(xsamp, 1, function(x) prox_nll(y, x, Bmat, x0))
  nll0 <- nll(y, x0, Bmat)
  mu <- as.numeric(-y * (t(Bmat) %*% x0))
  dn <- deta(mu); d2n <- d2eta(mu)
  gterm <- x0 + Bt %*% dn
  hterm <- eye(p) + Bt %*% diag(d2n) %*% t(Bt)
  delts <- t(t(xsamp) - x0)
  hdelts <- delts %*% hterm
  prox_nlls <- nll0 + as.numeric(delts %*% gterm) + 
    rowSums((delts %*% hterm) * delts)/2
  log_diffs <- prox_nlls - nlls
  (cf <- meanexp(log_diffs))
  se <- sd(exp(log_diffs))/sqrt(in.reps)
  if (log.p) {
    return(
      c(log.p = -l0 - log(det(hs))/2 + log(cf),
        cf = cf, se = se))
  }
  pr <- exp(-l0)/sqrt(det(hs)) * cf
  c(pr, se = se * pr)
}

compute_mi3 <- function(Bmat, mc.reps = 1e5, mcc = 0, in.reps = 1e2) {
  
  diffs <- unlist(mclapply0(1:mcc, function(i) {
    X <- randn(mc.reps, p)
    ps <- 1/(1 + exp(-X %*% Bmat))
    ces <- apply(ps, 1, function(v) sum(hp(v)))
    Y <- (rand(mc.reps, q) < ps) + 0
    logps <- apply(2 * Y - 1, 1, function(v) pr_laplace2(v, Bmat, TRUE, in.reps))
    -ces - logps[1, ]
  }, mc.cores = mcc))
  
  c(mu = mean(diffs), se = sd(diffs)/sqrt(length(diffs)))
}
