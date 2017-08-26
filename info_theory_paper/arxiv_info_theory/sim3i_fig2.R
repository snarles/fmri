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


mses <- colSums((res[, 1:7] - mi_true)^2)/data.reps

ihat <- as.numeric(res[, 1:6])
tab <- data.frame(I = ihat, method = rep(c("CM", "F", "LS", "0", "0.5", "0.9"), each = data.reps))
boxplot(I ~ method, data = tab)



boxplot(I ~ method, data = tab, ylim = c(0, 6),
        names = c(expression(hat(I)[0]), expression(hat(I)[0.5]), 
                  expression(hat(I)[0.9]), expression(hat(I)[CM]), 
                  expression(hat(I)[F]), expression(hat(I)[LS])),
        cex.axis = 1)
abline(mi_true, 0, lty = 2, col = "red", lwd = 3)
text(1.3, 6, "RMSE x 100 = ", col = "red", cex = 0.8)
reorder <- sqrt(c(mses[4:6], mses[1:3]))
for (i in 1:6) {
  text(i, 5.5, format(round(reorder[i] * 1e2, 1)),
       col = "red", cex = 0.8)
}
