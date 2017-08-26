source("info_theory_sims/sim3source.R")

####
##  Logistic simulation General B
####

allresults <- list()

## parallelization
mc.reps <- 1e5
mc.abe <- 1e4
mcc <- 39
data.reps <- mcc

## problem params
p <- 20; q <- 2 * p
Bmat <- 0.7/sqrt(p) * randn(p, q)
## true MI
# t1 <- proc.time()
# (mi_true <- compute_mi(Bmat, 1e6, mcc))
# proc.time() - t1
# 
# t1 <- proc.time()
# (mi_true2 <- compute_mi2(Bmat, 1e4, mcc))
# proc.time() - t1

t1 <- proc.time()
(mi_true3 <- compute_mi3(Bmat, 1e3, mcc, 1e3))
proc.time() - t1
#c(mi_true, mi_true2, mi_true3)
## bayes LS
k.each <- 20
(est_ls <- get_abe(Bmat, k.each, mc.abe, mcc))
## data params
m.folds <- 1
r.each <- 8000
r.train <- floor(0.5 * r.each)
(N = m.folds * k.each * r.each)
res <- run_simulations(Bmat, m.folds, k.each, r.each, r.train, mcc, data.reps)
## display results
c(mi_true = mi_true3, mi_ls = est_ls['mc_b_ls'], apply(res, 2, median), abe = est_ls['abe'])
apply(res[, 1:6] - mi_true, 2, summary)
colSums((res[, 1:6] - mi_true)^2)/data.reps
## save results
packet <- list(Bmat = Bmat, m.folds = m.folds,
               k.each = k.each, r.each = r.each, r.train = r.train,
               mi_true, est_ls = est_ls, res = res,
               mc.reps = mc.reps, mc.abe = mc.abe)
allresults <- c(allresults, list(packet))


save(allresults, file = 'info_theory_sims/fig3.Rdata')

####
##  plots
####
mi_true <- mi_true3[1]
mses <- colSums((res[, 1:7] - mi_true)^2)/data.reps

ihat <- as.numeric(res[, 1:6])
tab <- data.frame(I = ihat, method = rep(c("CM", "F", "LS", "0", "0.5", "0.9"), each = data.reps))

boxplot(I ~ method, data = tab, ylim = c(0, 5))

boxplot(I ~ method, data = tab, ylim = c(0, 5),
        names = c(expression(hat(I)[0]), expression(hat(I)[0.5]), 
                  expression(hat(I)[0.9]), expression(hat(I)[CM]), 
                  expression(hat(I)[F]), expression(hat(I)[LS])),
        cex.axis = 1)
abline(mi_true, 0, lty = 2, col = "red", lwd = 3)
text(1.3, 4.9, "RMSE x 100 = ", col = "red", cex = 0.8)
reorder <- sqrt(c(mses[4:6], mses[1:3]))
for (i in 1:6) {
  text(i, 4.35, format(round(reorder[i] * 1e2, 1)),
       col = "red", cex = 0.8)
}

