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
m.folds <- 1
r.each <- 8000
r.train <- floor(0.5 * r.each)
(N = m.folds * k.each * r.each)
t1 <- proc.time()
run_simulation(Bmat, m.folds, k.each, r.each, r.train)
proc.time() - t1


####
##   Generate an instance
####
set.seed(0)
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


nsubs = 10 * (1:800)
t1 <- proc.time()
res <- mclapply(nsubs, lala, mc.cores = mcc)
proc.time() - t1

res <- do.call(rbind, res)
head(res)
matplot(res[, 1:6], type = "l")
abline(mi_true, 0, lwd = 5)

res2 <- cbind(nsubs, res)
save(res2, file = "info_theory_sims/fig3data.RData")

sqrt(max(10 * nsubs))

plot(NA, NA, xlim = c(0, 290), ylim = c(0, 3), 
     xlab  = "n", ylab = "I", axes = FALSE)
axis(side = 1, labels = c("0", "10000", "40000","90000"),
     at = c(0, 100, 200, 300))
axis(side = 2)
lines(sqrt(10 * nsubs), res[, 'mi_0'], lwd = 6, col = gray(0.2))
lines(sqrt(10 * nsubs), res[, 'mi_0'], lwd = 3, col = "red")
lines(sqrt(10 * nsubs), res[, 'mi_5'], lwd = 6, col = gray(0.2))
lines(sqrt(10 * nsubs), res[, 'mi_5'], lwd = 3, col = "orange")
lines(sqrt(10 * nsubs), res[, 'mi_9'], lwd = 6, col = gray(0.2))
lines(sqrt(10 * nsubs), res[, 'mi_9'], lwd = 3, col = "yellow")
lines(sqrt(10 * nsubs), res[, 'mi_ls'], lwd = 6, col = gray(0.2))
lines(sqrt(10 * nsubs), res[, 'mi_ls'], lwd = 3, col = "green")
abline(mi_true, 0, col = "red", lty = 2, lwd = 3)
