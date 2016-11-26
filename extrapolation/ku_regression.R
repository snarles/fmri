library(lineId)

k <- 20
d <- 3 # degree

binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

get_avrisk_mat <- function(ks, d) {
  ans <- matrix(0, length(ks), d + 1)
  ans <- col(ans) - 1
  ans <- ans + (ks - 1)
  ans <- 1/ans
  ans <- ans * (ks - 1)
  ans
}

get_kmat <- function(k, d) {
  kmat <- matrix(0, k - 1, d + 1)
  al <- row(kmat)
  bt <- k - al
  h <- col(kmat) - 1
  kmat <- lgamma(al + h) + lgamma(al + bt)  - lgamma(al) - lgamma(al + bt + h)
  kmat <- exp(kmat)
  kmat  
}

get_rank_prop <- function(pmat, true_ys) {
  k <- nrow(pmat)
  p2 <- apply(pmat, 2, rank)
  true_ranks <- p2[cbind(true_ys, 1:ncol(pmat))]
  tab <- table(true_ranks)
  ans <- numeric(k)
  ans[as.numeric(names(tab))] <- tab
  ans <- ans/ncol(pmat)
  cumsum(ans)[1:(k-1)]
}

get_sub_errs <- function(pmat, true_ys, ks) {
  k <- nrow(pmat)
  p2 <- apply(pmat, 2, rank)
  true_ranks <- p2[cbind(true_ys, 1:ncol(pmat))]
  1 - sapply(ks, function(v) mean(binmom(true_ranks - 1, k - 1, v - 1)))
}

get_vande <- function(xs = seq(0, 1, 0.01), d) {
  ans <- xs %*% t(rep(1, d+1)) 
  ans <- ans ^ (col(ans) - 1)
  ans
}

avrisk <- function(k, bt) {
  d <- length(bt) - 1
  sapply(k, function(k) (k - 1) * sum(bt/(0:d + k - 1)))
}

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
