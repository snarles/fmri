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
  st/sds
}

double_boot <- function(boot_stat, XA, YA, XB, YB, rep.out = 1000, rep.in = rep.out) {
  nA <- dim(YA)[1]; nB <- dim(YB)[1]
  outer.vals <- numeric()
  for (i in 1:rep.out) {
    indA <- sample(nA, nA, replace = TRUE)
    indB <- sample(nB, nB, replace = TRUE)
    outer.vals <- rbind(outer.vals,
                        inner_boot(boot_stat, XA[indA, ], YA[indA, ], XB[indB, ], YB[indB, ], rep.in))
  }
  outer.vals
}

####
## TEST
####

nA <- 40; nB <- 50; p <- 2; q <- 10
XA <- randn(nA, p)
XB <- randn(nB, p)
G0 <- randn(p, q); G1 <- randn(p, q)
GA <- svd(randn(q, q))$u
GB <- svd(randn(q, q))$u
# null case
BA0 <- G0 %*% GA; BB0 <- G0 %*% GB
# non-null case
BA0 <- G0 %*% GA; BB0 <- G1 %*% GB

SigmaA <- cov(randn(2 * q, q))
SigmaB <- cov(randn(2 * q, q))
truth <- list(G=G0 %*% t(G0), BA = BA0, BB=BB0)


YA <- XA %*% BA0 + mvrnorm(nA, rep(0, q), SigmaA)
YB <- XB %*% BB0 + mvrnorm(nB, rep(0, q), SigmaA)

boot_stat(XA, YA, XB, YB)

inner_boot(boot_stat, XA, YA, XB, YB)
vals <- double_boot(boot_stat, XA, YA, XB, YB, rep.out = 100)
layout(matrix(1:4, 2, 2))
hist(vals[, 1])
hist(vals[, 2])
hist(vals[, 3])
hist(vals[, 4])

