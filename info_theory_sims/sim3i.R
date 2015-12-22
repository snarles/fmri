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
p <- 50
mult <- 4/sqrt(p)
Bmat <- mult * eye(p)
(mi_true <- mi_ident_case(p, mult, 1e5))
## bayes LS
k.each <- 10
t1 <- proc.time()
(est_ls <- get_abe(Bmat, k.each, abe.each, mc.abe, mcc))
proc.time() - t1
## data params
m.folds <- 1
r.each <- 8000
r.train <- floor(0.5 * r.each)
(N = m.folds * k.each * r.each)
t1 <- proc.time()
run_simulation(Bmat, m.folds, k.each, r.each, r.train)
proc.time() - t1
## full-scale
res <- run_simulations(Bmat, m.folds, k.each, r.each, r.train, mcc, data.reps)
## display results
c(mi_true = mi_true, mi_ls = est_ls['mc_b_ls'], apply(res, 2, median), abe = est_ls['abe'])
apply(res[, 1:6] - mi_true, 2, summary)
colSums((res[, 1:6] - mi_true)^2)/data.reps
## save results
packet <- list(Bmat = Bmat, m.folds = m.folds,
               k.each = k.each, r.each = r.each, r.train = r.train,
               mi_true, est_ls = est_ls, res = res,
               mc.reps = mc.reps, mc.abe = mc.abe)
allresults <- c(allresults, list(packet))


save(allresults, file = 'info_theory_sims/save.Rdata')