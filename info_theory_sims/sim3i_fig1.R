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
data.reps <- 90

## problem params
p <- 3
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
r.each <- 40
r.train <- floor(0.5 * r.each)
(N = m.folds * k.each * r.each)
t1 <- proc.time()
run_simulation(Bmat, m.folds, k.each, r.each, r.train)
proc.time() - t1
## full-scale
t1 <- proc.time()
res <- run_simulations(Bmat, m.folds, k.each, r.each, r.train, mcc, data.reps)
proc.time() - t1
## display results
c(mi_true = mi_true, mi_ls = est_ls['mc_b_ls'], apply(res, 2, median), abe = est_ls['abe'])
apply(res[, 1:6] - mi_true, 2, summary)
mses <- colSums((res[, 1:7] - mi_true)^2)/data.reps

ihat <- as.numeric(res[, 1:6])
tab <- data.frame(I = ihat, method = rep(c("CM", "F", "LS", "0", "0.5", "0.9"), data.reps))
View(tab)

boxplot(I ~ method, data = tab, ylim = c(0, 1.5))

boxplot(I ~ method, data = tab, ylim = c(0, 1.5),
        names = c(expression(hat(I)[0]), expression(hat(I)[0.5]), 
                  expression(hat(I)[0.9]), expression(hat(I)[CM]), 
                  expression(hat(I)[F]), expression(hat(I)[LS])),
        cex.axis = 1)
abline(mi_true, 0, lty = 2, col = "red", lwd = 3)
text(1.3, 1.5, "RMSE x 100 = ", col = "red", cex = 0.8)
reorder <- sqrt(c(mses[4:6], mses[1:3]))
for (i in 1:6) {
  text(i, 1.35, format(round(reorder[i] * 1e2, 1)),
       col = "red", cex = 0.8)
}
