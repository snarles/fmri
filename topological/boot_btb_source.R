library(pracma)
library(MASS)
library(lineId)
library(parallel)
#library(ade4)

####
##  Double bootstrap method
####

boot_stat <- function(XA, YA, XB, YB) {
  BA <- solve(t(XA) %*% XA, t(XA) %*% YA)
  BB <- solve(t(XB) %*% XB, t(XB) %*% YB)
  as.numeric(BA %*% t(BA) - BB %*% t(BB))
}

inner_boot <- function(boot_stat, XA, YA, XB, YB, rep.in = 1000) {
  vals <- numeric()
  st <- boot_stat(XA, YA, XB, YB)
  for (i in 1:rep.in) {
    indA <- sample(nA, nA, replace = TRUE)
    indB <- sample(nB, nB, replace = TRUE)
    vals <- rbind(vals, boot_stat(XA[indA, ], YA[indA, ], XB[indB, ], YB[indB, ]))
  }
  sds <- apply(vals, 2, sd)
  ## Without bias-correction
  #st/sds
  ## With bias-correction
  (2 * st - colMeans(vals))/sds
}

double_boot <- function(boot_stat, XA, YA, XB, YB, rep.out = 1000,
                        rep.in = rep.out, mcc = 0) {
  nA <- dim(YA)[1]; nB <- dim(YB)[1]
  lala <- function (i) {
    indA <- sample(nA, nA, replace = TRUE)
    indB <- sample(nB, nB, replace = TRUE)
    inner_boot(boot_stat, XA[indA, ], YA[indA, ], XB[indB, ], YB[indB, ], rep.in)
  }  
  outer.vals <- mclapply0(1:rep.out, lala, mc.cores=mcc)
  do.call(rbind, outer.vals)
}

veczero_pvalue <- function(vals) {
  p <- dim(vals)[2]
  n <- dim(vals)[1]
  qt <- apply(vals, 2, function(v) sum(v > 0)/length(v))
  pvals <- 2 * pmin(qt, 1-qt)
  bonf_p <- min(pvals) * p
  bonf_p
}

####
## TEST
####

# nA <- 40; nB <- 50; p <- 2; q <- 10
# XA <- randn(nA, p)
# XB <- randn(nB, p)
# G0 <- randn(p, q); G1 <- randn(p, q)
# GA <- svd(randn(q, q))$u
# GB <- svd(randn(q, q))$u
# # null case
# BA0 <- G0 %*% GA; BB0 <- G0 %*% GB
# # non-null case
# BA0 <- G0 %*% GA; BB0 <- G1 %*% GB
# 
# SigmaA <- cov(randn(2 * q, q))
# SigmaB <- cov(randn(2 * q, q))
# truth <- list(G=G0 %*% t(G0), BA = BA0, BB=BB0)
# 
# 
# YA <- XA %*% BA0 + mvrnorm(nA, rep(0, q), SigmaA)
# YB <- XB %*% BB0 + mvrnorm(nB, rep(0, q), SigmaA)
# 
# boot_stat(XA, YA, XB, YB)
# 
# inner_boot(boot_stat, XA, YA, XB, YB)
# vals <- double_boot(boot_stat, XA, YA, XB, YB, rep.out = 100)
# layout(matrix(1:4, 2, 2))
# hist(vals[, 1])
# hist(vals[, 2])
# hist(vals[, 3])
# hist(vals[, 4])
# 
# veczero_pvalue(vals)
