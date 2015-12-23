source("info_theory_sims/sim3source.R")

####
##  Identity case
####

allresults <- list()

## parallelization
mc.reps <- 1e4
mc.abe <- 1e2
abe.each <- 1e2
mcc <- 3
data.reps <- mcc

## problem params
p <- 10
mult <- 2*sqrt(2)/sqrt(p)
Bmat <- mult * eye(p)
(mi_true <- mi_ident_case(p, mult, 1e5))
## bayes LS
k.each <- 5
t1 <- proc.time()
(est_ls <- get_abe(Bmat, k.each, abe.each, mc.abe, mcc))
proc.time() - t1
## data params
m.folds <- 1
r.each <- 40000
r.train <- floor(0.5 * r.each)
(N = m.folds * k.each * r.each)
t1 <- proc.time()
run_simulation(Bmat, m.folds, k.each, r.each, r.train)
proc.time() - t1


####
##   Generate an instance
####
q <- dim(Bmat)[2]
X <- randn(k.each, p)
ps <- 1/(1 + exp(-X %*% Bmat))
zs <- rep(1:k.each, r.each)
Yall <- (rand(k.each * r.each, q) < ps[zs, ]) + 0

nsub <- k.each * 10
ntr <- k.each * 5

lala <- function(i) {
  run_instance(X, Yall, zs, 10 * i, 5 * i)
}

lala(3)
