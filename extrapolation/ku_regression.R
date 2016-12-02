source("extrapolation/ku_source.R")

## get pmat
load("rakesh/converted1.rda", verbose = TRUE)
# err1$TestErr <- 1 - err1$TestErr
#View(err1)
#View(err_knn400)
lprobs2 <- readRDS("rakesh/converted2.rds")


(ac0 <- 1 - err1[err1$configuration=="logistic" & err1$num_classes==400, "TestErr"])
pmat <- lprobs$logistic_111938.logprobs
dim(pmat) # 20 1000
true_ys <- rep(1:20, each = 50)

####
##  polynomial model
####


k <- nrow(pmat)
d <- 5

ws <- exp(0 * (1:(k-1)))
ws <- ws/sum(ws)
kmat <- get_kmat(k, d)
cvec <- get_rank_prop(pmat, true_ys)

(bt <- solve(t(kmat) %*% diag(ws) %*% kmat, t(kmat) %*% diag(ws) %*% cvec))
sum(bt)

ku <- get_vande(d = d) %*% bt


yhat <- kmat %*% bt
yhat

plot(kmat[, 2], cvec, ylim = c(0, 0.5))
lines(kmat[, 2], yhat, col = "red")
lines(seq(0, 1, 0.01), ku, col = "blue")

avrisk(k, bt)
cvec
avrisk(30, bt)
plot(2:200, avrisk(2:200, bt), type = "l")


####
##  Constraining predictions
####

library(limSolve)
k <- nrow(pmat)
d <- 5
kmat <- get_kmat(k, d)
cvec <- get_rank_prop(pmat, true_ys)

#ks <- c(k - 2, k -1 , k)
ks <- (k-5):k
avr_w <- get_avrisk_mat(ks, d)
avr <- get_sub_errs(pmat, true_ys, ks)

bt <- lsei(A = kmat, B = cvec)$X
sum(bt)

bt <- lsei(A = kmat, B = cvec, E = avr_w, F = avr)$X
sum(bt)

Ws <- exp(5 * (2:k)/k)
bt <- lsei(A = Ws * get_avrisk_mat(2:k, d), B = Ws * get_sub_errs(pmat, true_ys, 2:k))$X
sum(bt)

ku <- get_vande(d = d) %*% bt
yhat <- kmat %*% bt
plot(kmat[, 2], cvec, ylim = c(0, 0.5))
lines(kmat[, 2], yhat, col = "red")
lines(seq(0, 1, 0.01), ku, col = "blue")

plot(2:200, avrisk(2:200, bt), type = "l")


####
##  Spline basis
####

nsplines <- 20
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
ffs <- lapply(knts, spline1_maker)

Kmax <- 200

ks <- (k-5):k
MM <- make_moment_mat(ffs, 1:Kmax, res = 1e5)
xmat <- MM[ks, ]
avr <- get_sub_errs(pmat, true_ys, ks)

#bt <- pinv(xmat) %*% avr
bt <- nnls(xmat, avr)$X
plot(MM %*% bt, type = "l")
1-ac0
