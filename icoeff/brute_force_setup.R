## gets matrix pairs plus calculates power spectra up to order deg
library(parallel)
library(abind)
source("icoeff/reg_graphon.R")

set.seed(0)
mc.cores <- 40
p <- 5
self.loop <- TRUE
normalize <- TRUE
shape <- 1
rate <- 1
n.pairs <- 10
eps <- 0.1 ## make the perturbation small
deg <- 12

before_after_pair <- function(sd) {
  set.seed(sd)
  P <- gen_rwg(shape, rate, p = p, self.loop = self.loop, normalize = normalize)
  S <- (1 - eps) * eye(p) + eps * gen_rwg(shape, rate, p = p, self.loop = self.loop, normalize = FALSE)
  #check_reg(S)
  rbind(P, P %*% S, S %*% P)
}

mi <- function(mat) {
  #stopifnot(sum(mat) == 1)
  sum(mat *log(mat)) + 2 *log(nrow(mat))
}
before_after_pair_mi <- function(sd) {
  set.seed(sd)
  P <- gen_rwg(shape, rate, p = p, self.loop = self.loop, normalize = normalize)
  S <- (1 - eps) * eye(p) + eps * gen_rwg(shape, rate, p = p, self.loop = self.loop, normalize = FALSE)
  #check_reg(S)
  c(mi(P), mi(P %*% S), mi(S %*% P))
}


matblock <- mclapply(1:n.pairs, before_after_pair, mc.cores = mcc)
mis <- mclapply(1:n.pairs, before_after_pair_mi, mc.cores = mcc)
res_mi <- do.call(cbind, mis)

## iterated powering

m0 <- do.call(rbind, matblock)
dim(m0)
m <- m0
sumsS <- list()

for (i in 1:(deg) + 1) {
  raws <- rowSums(m)
  sums <- colSums(matrix(raws, p, n.pairs * 3))
  sums <- matrix(sums, 3, n.pairs)
  sumsS[[i]] <- sums
  m <- m * m0
}

help(abind)
a_sums <- abind(sumsS, along = 3)
a_sums <- aperm(a_sums, c(3, 1, 2))

## this will be moved to another script
coefs <- readRDS("icoeff/coeffs_xlogx.rds")
coefs <- coefs[1:(deg + 1)]
evald <- coefs[1] + array_colsum(a_sums * coefs[-1]) + 2 * log(p)
dim(evald)
evald



mi2 <- function(mat) {
  #stopifnot(sum(mat) == 1)
  m <- mat
  s <- coefs[1]
  for (i in 2:(deg + 1)) {
    s <- s + coefs[i] * sum(m)
    m <- m * mat
  }
  s + 2 * log(nrow(mat))
}
P <- .9 * eye(5) + .1 * gen_rwg(p = 5, normalize = TRUE)
mi(P)
mi2(P)
