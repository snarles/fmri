library(lineId)

k <- 20
d <- 3 # degree

get_kmat <- function(k, d) {
  kmat <- matrix(0, k - 1, d + 1)
  current <- rep(1, k - 1)
  for (i in 0:d) {
    kmat[, i + 1] <- current
    current <- current * ((0 : (k-2)) - i)/(k - 2 - i)
  }
  kmat[, i + 1] <- current
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

k <- nrow(pmat)
d <- 5

ws <- exp(0 * (1:(k-1)))
ws <- ws/sum(ws)
kmat <- get_kmat(k, d)
cvec <- get_rank_prop(pmat, true_ys)

(bt <- solve(t(kmat) %*% diag(ws) %*% kmat, t(kmat) %*% diag(ws) %*% cvec))

ku <- get_vande(d = d) %*% bt

yhat <- kmat %*% bt

plot(kmat[, 2], cvec, ylim = c(0, 0.5))
lines(kmat[, 2], yhat, col = "red")
lines(seq(0, 1, 0.01), ku, col = "blue")

avrisk(k, bt)
cvec
avrisk(30, bt)
plot(2:200, avrisk(2:200, bt), type = "l")
