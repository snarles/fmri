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
data.reps <- 39

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
t1 <- proc.time()
#run_simulation(Bmat, m.folds, k.each, r.each, r.train)
proc.time() - t1
## full-scale
t1 <- proc.time()
#res <- run_simulations(Bmat, m.folds, k.each, r.each, r.train, mcc, data.reps)
load("info_theory_sims/fig2.Rdata")
proc.time() - t1

## compute abe computation

library(lineId)
ident <- sapply(1 - res[, "abe"], function(v) aba_to_mi_lower(k.each, v))


## display results
c(mi_true = mi_true, mi_ls = est_ls['mc_b_ls'], apply(res, 2, median), abe = est_ls['abe'])
apply(res[, 1:6] - mi_true, 2, summary)
mses <- colSums((res[, 1:7] - mi_true)^2)/data.reps
mse_ident <- sum((ident - mi_true)^2)/data.reps

ihat <- c(as.numeric(res[, 1:6]), ident)
tab <- data.frame(I = ihat, method = rep(c("CM", "F", "LS", "0", "0.5", "0.9", "Id"), each = data.reps))

boxplot(I ~ method, data = tab, ylim = c(0, 6))

boxplot(I ~ method, data = tab, ylim = c(0, 6),
        names = c(expression(hat(I)[0]), expression(hat(I)[0.5]), 
                  expression(hat(I)[0.9]), expression(hat(I)[CM]), 
                  expression(hat(I)[F]), expression(hat(I)[Id]),
                  expression(hat(I)[HD])),
        cex.axis = 1)
abline(mi_true, 0, lty = 2, col = "red", lwd = 3)
text(1.5, 6, "RMSE x 100 = ", col = "red", cex = 0.8)
reorder <- sqrt(c(mses[4:6], mses[1:2], mse_ident, mses[3]))
for (i in 1:7) {
  text(i, 5.5, format(round(reorder[i] * 1e2, 1)),
       col = "red", cex = 0.8)
}
#save(res, file = "info_theory_sims/fig1.Rdata")
