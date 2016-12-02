source("extrapolation/ku_source.R")

k <- 20
d <- 3 # degree





true_d <- 3
k <- 300
sigma <- 1
mu <- randn(k, true_d)
y <- mu + sigma * randn(k, true_d)
pmat <- -pdist2(mu, y)
true_ys <- 1:k

####
##  polynomial model
####


k <- nrow(pmat)
d <- 10

ws <- exp(0 * (1:(k-1)))
ws <- ws/sum(ws)
kmat <- get_kmat(k, d)
cvec <- get_rank_prop(pmat, true_ys)

(bt <- solve(t(kmat) %*% kmat, t(kmat) %*% cvec))
(bt <- solve(t(kmat) %*% diag(ws) %*% kmat, t(kmat) %*% diag(ws) %*% cvec))
sum(bt)

ku <- get_vande(d = d) %*% bt


yhat <- kmat %*% bt
yhat

plot(kmat[, 2], cvec, ylim = c(0, 1), type = "l")
lines(kmat[, 2], yhat, col = "red")
lines(seq(0, 1, 0.01), ku, col = "blue")

avrisk(k, bt)
cvec
avrisk(5, bt)

avrisk(30, bt)
avrisk(100, bt)
plot(20:300, avrisk(20:300, bt), type = "l")


####
##  Constraining predictions
####

library(limSolve)
k <- nrow(pmat)
d <- 10
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

Ws <- exp(0 * (2:k)/k)
bt <- lsei(A = Ws * get_avrisk_mat(2:k, d), B = Ws * get_sub_errs(pmat, true_ys, 2:k))$X


ks <- 20:k
avr_w <- get_avrisk_mat(ks, d)
avr <- get_sub_errs(pmat, true_ys, ks)
(bt <- pinv(avr_w) %*% avr)
sum(bt)

ku <- get_vande(d = d) %*% bt
yhat <- kmat %*% bt
plot(kmat[, 2], cvec, ylim = c(0, 1), type = "l")
lines(kmat[, 2], yhat, col = "red")
lines(seq(0, 1, 0.01), ku, col = "blue", type = "l")

plot(20:300, avrisk(20:300, bt), type = "l")


#####
##  Spline
#####

nsplines <- 20
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
ffs <- lapply(knts, spline1_maker)

Kmax <- 500

ks <- 20:k
MM <- make_moment_mat(ffs, 1:Kmax, res = 1e5)
xmat <- MM[ks, ]
avr <- get_sub_errs(pmat, true_ys, ks)

#bt <- pinv(xmat) %*% avr
bt <- nnls(xmat, avr)$X
sum(bt)
plot(MM %*% bt, type = "l")

ts <- seq(0, 1, 0.01)
xs <- sapply(ffs, function(ff) ff(ts))
ku <- as.numeric(xs %*% bt)

plot(kmat[, 2], cvec, ylim = c(0, 1), type = "l")
lines(seq(0, 1, 0.01), ku, col = "blue", type = "l")




