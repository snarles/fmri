library(pracma)
library(MASS)
library(class)
library(parallel)
library(reginference) ## see github.com/snarles/misc
library(lineId) ## use devtools::install('lineId')
#source("approximation/extreme_value.R")


## log 1+x
log_1_plus <- function(x) {
  if (abs(x) > 1e-3) return(log(1 + x))
  sum(-1 * (-x)^(1:20)/(1:20))
}

## mu ~ N(0, I)
## y ~ N(mu*, sigma^2 I)

## most naive implementation of mc() for identity cov
mc_ident <- function(p, sigma2, K, mc.reps = 1000) {
  mcs <- sapply(1:mc.reps,
                function(i) {
                  mus <- randn(K, p)
                  ys <- mus + sqrt(sigma2) * randn(K, p)
                  lbls <- knn(mus, ys, cl=1:K)
                  sum(lbls != 1:K)/K
                })
  mean(mcs)
}

## uses noncentral chi squared
## same Evalue as mc_ident
mc_ident2 <- function(p, sigma2, K, mc.reps = 1000) {
  alpha <- sigma2/(1 + sigma2)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y2 <- (1 + sigma2) * rchisq(1, df = p)
                  d1 <- alpha * rchisq(1, df = p, ncp = alpha * y2)
                  ds <- rchisq(K - 1, df = p, ncp = y2)
                  (min(ds) < d1)
                })
  mean(mcs)
}



## actually generates a gaussian instead of a chi-squared
rchisq_g <- function(n, df, ncp = 0) {
  mu <- df + ncp
  vv <- 2 * df + 4 * ncp
  sqrt(vv) * rnorm(n) + mu
}

## uses gaussian--approximate!!
mc_ident3 <- function(p, sigma2, K, mc.reps = 1000) {
  alpha <- sigma2/(1 + sigma2)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y2 <- (1 + sigma2) * rchisq_g(1, df = p)
                  if (y2 > 0) {
                    d1 <- alpha * rchisq_g(1, df = p, ncp = alpha * y2)
                    ds <- rchisq_g(K - 1, df = p, ncp = y2)
                    return(min(ds) < d1)
                  } else {
                    return(1 - (1/K)) ## default in case of error
                  }
                })
  mean(mcs)
}

min_chisq_approx <- function(K, df, ncp, naive = FALSE, nits = 20) {
  if (naive) {
    print(min(rchisq(K, df, ncp)))
  }
  #z <- (ncp - df)/sqrt(df)
  #med_loc2(df, z, K, ux=df + ncp, nits=nits)
  ux <- df + ncp; lx <- 0
  const <- -log(2)
  for (i in 1:nits) {
    xc <- (ux + lx)/2
    val <- K * log_1_plus(-pchisq(xc, df, ncp))
    if (val > const) {
      lx <- xc
    } else {
      ux <- xc
    }
    #print(xc)
  }
  xc
}

## computes theoretical median of min
mc_ident4 <- function(p, sigma2, K, mc.reps = 1000, nits = 20) {
  alpha <- sigma2/(1 + sigma2)
  mcs <- sapply(1:mc.reps,
                function(i) {
                  y2 <- (1 + sigma2) * rchisq(1, df = p)
                  d1 <- alpha * rchisq(1, df = p, ncp = alpha * y2)
                  ds <- min_chisq_approx(K - 1, p, y2, nits = nits)
                  (ds < d1)
                })
  mean(mcs)
}

## nonrandom -- computes theoretical median of min and uses cdf
## same Evalue as mc_ident4
mc_ident5 <- function(p, sigma2, K, mc.reps = 1000, nits = 20) {
  alpha <- sigma2/(1 + sigma2)
  quants <- (1:mc.reps - .5)/mc.reps
  y2s <- (1 + sigma2) * qchisq(quants, df = p)
  mcs <- sapply(y2s,
                function(y2) {
                  ds <- min_chisq_approx(K - 1, p, y2, nits = nits)
                  1 - pchisq(ds/alpha, p, alpha * y2) 
                })
  mean(mcs)
}



# mean(rchisq(1e3, 3, 2.2))
# mean(rchisq_g(1e3, 3, 2.2))
# var(rchisq(1e3, 3, 2.2))
# var(rchisq_g(1e3, 3, 2.2))
# 
# ## Gaussian approximation is unusably bad
# 
# mc_ident(10, 3, 10, 1e4)
# mc_ident2(10, 3, 10, 1e4)
# mc_ident4(10, 3, 10, 1e4)
# 
# mc_ident( 20, 1, 1e3, 1e3)
# mc_ident2(20, 1, 1e3, 1e3)
# mc_ident4(20, 1, 1e3, 1e3)
# 
# 
# df <- 40; sigma2 <- 2; L <- 1e5; nits <- 100
# 
# tims <- list(); res <- numeric()
# t1 <- proc.time()
# (res[1] <- mc_ident( df, sigma2, L, 1))
# (tims[[1]] <- proc.time() - t1)
# 
# t1 <- proc.time()
# (res[2] <- mc_ident2(df, sigma2, L, nits))
# (tims[[2]] <- proc.time() - t1)
# 
# t1 <- proc.time()
# (res[3] <- mc_ident4(df, sigma2, L, nits))
# (tims[[3]] <- proc.time() - t1)
# 
# res
# tims
# 
# 
# ## push the dimensionality and # of classes
# 
# min_chisq_approx(1e5, 20, 20, TRUE)
# min_chisq_approx(1e6, 30, 30, TRUE)
# t1 <- proc.time()
# min_chisq_approx(1e8, 40, 40, TRUE); gc()
# proc.time() - t1


