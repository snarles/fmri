#####
##  Checking out model: monotonic g model
#####

k <- 7
dgu <- runif(k + 1)
gu <- cumsum(dgu)
(gu <- gu/sum(gu))
(dgu <- c(gu[1], gu[-1] - gu[-(k + 1)]))
# rbind(cumsum(dgu), gu)
binprobs <- matrix(0, k + 1, k + 1)
for (i in 0:k) binprobs[, i + 1] <- dbinom(0:k, k, i/k)
# matplot(binprobs, type = "l")
hy <- binprobs %*% gu
vprobs <- matrix(0, k + 1, k + 1)
for (i in 0:k) vprobs[, i+1] <- binprobs %*% c(rep(0, i), rep(1, k + 1 - i))
hy2 <- vprobs %*% dgu
cbind(hy, hy2)

lambda <- 2
bt <- runif(k + 1); bt <- bt/sum(bt)
(bt <- bt[-(k+1)])
bt0 <- bt
ft <- vprobs %*% c(bt, 1-sum(bt))
(ys <- sample(k + 1, 20, TRUE, prob = ft))
(ws <- sapply(0:k, function(i) sum(ys == i)))
bt <- runif(k + 1); bt <- bt/sum(bt)
(bt <- bt[-(k+1)])
of <- function(bt) {
  ft <- vprobs %*% c(1-sum(bt), bt)
  -sum(ws * log(ft)) - lambda * (sum(log(bt)) + log(1 - sum(bt)))
}
og <- function(bt) {
  k <- length(bt)
  ans <- numeric(k)
  ft <- vprobs %*% c(1-sum(bt), bt)
  for (i in 1:k) {
    s <- -sum(ws * (vprobs[, i + 1] - vprobs[, 1])/ft)
    s <- s - lambda/bt[i] + lambda/(1-sum(bt))
    ans[i] <- s
  }
  ans
}
oh <- function(bt) {
  k <- length(bt)
  ans <- matrix(0, k, k)
  ft <- vprobs %*% c(1-sum(bt), bt)
  for (i in 1:k) {
    for (j in 1:k) {
      s <- sum(ws * (vprobs[, i + 1] - vprobs[, 1]) * 
                 (vprobs[, j + 1] - vprobs[, 1])/(ft^2))
      if (i == j) s <- s + lambda/(bt[i]^2)
      s <- s + lambda/(1 - sum(bt))^2
      ans[i, j] <- s
    }
  }
  ans
}
of(bt)
of(bt0)
numDeriv::grad(of, bt0)
og(bt0)
numDeriv::hessian(of, bt0)
oh(bt0)
