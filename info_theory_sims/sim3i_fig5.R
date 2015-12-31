source("info_theory_sims/sim3source.R")

####
##  Identity case
####

## parallelization
mc.reps <- 1e4
mc.abe <- 1e2
abe.each <- 1e2
mcc <- 39
data.reps <- mcc * 3

## problem params
p <- 10
mult <- 4/sqrt(p)
Bmat <- mult * eye(p)
(mi_true <- mi_ident_case(p, mult, 1e5))
## bayes LS
k.each <- 20
t1 <- proc.time()
(est_ls <- get_abe(Bmat, k.each, abe.each, mc.abe, mcc))
proc.time() - t1
## data params
N <- 160000
ks <- c(5, 10, 15, 20, 25, 30, 35, 40)
ress <- lapply(ks, function(k) 
  run_simulations(Bmat, 1, k, floor(N/k), floor(0.5 * N/k), mcc, data.reps))
#res <- run_simulations(Bmat, m.folds, k.each, r.each, r.train, mcc, data.reps)

save(ress, file = "info_theory_sims/fig5.Rdata")

res2 <- do.call(rbind, ress)
K <- rep(ks, each = data.reps)
tab1 <- data.frame(I = res2[, 3], K = K)
res3 <- do.call(cbind, ress)
mses <- colSums((res3 - mi_true)^2)/data.reps
mses_ls <- mses[(1:length(ks)) * 9 - 6]

boxplot(I ~ K, data = tab1, ylim = c(0.5, 2.8))
abline(mi_true, 0, lty = 2, col = "red", lwd = 3)
text(1, 2.8, "RMSE x 100 = ", col = "red", cex = 0.8)
for (i in 1:length(ks)) text(i, 2.5, format(round(mses_ls[i] * 1e2, 1)), col = "red")
