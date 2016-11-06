source("info_theory_sims/sim4source.R")

####
##  Identity case
####

allresults <- list()

## parallelization
mc.reps <- 1e4
mc.abe <- 1e2
abe.each <- 1e2
mcc <- 3
data.reps <- 12

## problem params
p <- 8
mult <- 10/sqrt(p)
Bmat <- mult * eye(p)
(mi_true <- mi_ident_case(p, mult, 1e5))
## bayes LS
k.each <- 40
t1 <- proc.time()
(est_ls <- get_abe(Bmat, k.each, abe.each, mc.abe, mcc))
proc.time() - t1
## data params
m.folds <- 5
r.each <- 200
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

ihat <- as.numeric(res[, c(1:6, 11, 12)])
tab <- data.frame(I = ihat, 
                  method = rep(c("CM", "F", "LS", "0", "0.5", "0.9", "LB", "L2"), 
                               each = data.reps))

boxplot(I ~ method, data = tab)

# boxplot(I ~ method, data = tab, ylim = c(0, 1.5),
#         names = c(expression(hat(I)[0]), expression(hat(I)[0.5]), 
#                   expression(hat(I)[0.9]), expression(hat(I)[CM]), 
#                   expression(hat(I)[F]), expression(hat(I)[LS])),
#         cex.axis = 1)
abline(mi_true, 0, lty = 2, col = "red", lwd = 3)
# text(1.3, 3, "RMSE x 100 = ", col = "red", cex = 0.8)
# reorder <- sqrt(c(mses[4:6], mses[1:3]))
# for (i in 1:6) {
#   text(i, 2.8, format(round(reorder[i] * 1e2, 1)),
#        col = "red", cex = 0.8)
# }
save(res, file = "info_theory_sims/sim4fig1.Rdata")
