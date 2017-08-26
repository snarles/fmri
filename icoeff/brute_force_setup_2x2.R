## gets matrix pairs plus calculates power spectra up to order deg
## CONJECTURE: for p = 2, all cycles decrease

library(parallel)
library(abind)
source("icoeff/reg_graphon.R")

set.seed(0)
mc.cores <- 40
p <- 2
self.loop <- TRUE
normalize <- TRUE
shape <- 1
rate <- 1
n.pairs <- 1e4
eps <- 1e-3 ## make the perturbation small
deg <- 12

before_after_pair <- function(sd) {
  set.seed(sd)
  P <- gen_rwg(shape, rate, p = p, self.loop = self.loop, normalize = normalize)
  S <- (1 - eps) * eye(p) + eps * gen_rwg(shape, rate, p = p, self.loop = self.loop, normalize = FALSE)
  #check_reg(S)
  #rbind(P, P %*% S, S %*% P)
  rbind(P, P %*% S)
}

t1 <- proc.time()
matblock <- mclapply(1:n.pairs, before_after_pair, mc.cores = mcc)
proc.time() - t1

## iterated powering

m0 <- do.call(rbind, matblock)
dim(m0)
m <- m0
sumsS <- list()

t2 <- proc.time()
for (i in 1:(deg) + 1) {
  raws <- rowSums(m)
  sums <- colSums(matrix(raws, p, n.pairs * 2))
  sums <- matrix(sums, 2, n.pairs)
  sumsS[[i]] <- sums
  m <- m * m0
}
proc.time() - t2

a_sums <- abind(sumsS, along = 3)
a_sums <- aperm(a_sums, c(3, 1, 2))

a_diffs <- a_sums[,1,] - a_sums[,2,]
a_diffs <- a_diffs[-1, ]

sum(a_diffs < 0)
