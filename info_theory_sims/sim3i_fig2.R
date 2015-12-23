source("info_theory_sims/sim3source.R")

####
##  Identity case
####

## parallelization
mc.reps <- 1e4
mc.abe <- 1e2
abe.each <- 1e2
mcc <- 39
data.reps <- mcc

## problem params
p <- 50
mult <- 4/sqrt(p)
Bmat <- mult * eye(p)
(mi_true <- mi_ident_case(p, mult, 1e5))
## bayes LS
k.each <- 20
t1 <- proc.time()
(est_ls <- get_abe(Bmat, k.each, abe.each, mc.abe, mcc))
proc.time() - t1
## data params
m.folds <- 1
r.each <- 8000
r.train <- floor(0.5 * r.each)
(N = m.folds * k.each * r.each)
# t1 <- proc.time()
# run_simulation(Bmat, m.folds, k.each, r.each, r.train)
# proc.time() - t1
## full-scale
res <- run_simulations(Bmat, m.folds, k.each, r.each, r.train, mcc, data.reps)
save(res, file = "info_theory_sims/fig2.Rdata")
