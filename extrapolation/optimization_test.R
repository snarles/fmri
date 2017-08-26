####
##  Test numerical stability of stuff
####

l1d <- function(x, y = 0) sum(abs(x - y))

make_vprob <- function(k) {
  binprobs <- matrix(0, k + 1, k + 1)
  for (i in 0:k) binprobs[, i + 1] <- dbinom(0:k, k, i/k)
  # matplot(binprobs, type = "l")
  vprobs <- matrix(0, k + 1, k + 1)
  for (i in 0:k) vprobs[, i+1] <- (binprobs %*% c(rep(0, i), rep(1, k + 1 - i)))/(k + 1 - i)
  vprobs
}

make_functions <- function(vprobs, ws, lambda) {
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
  list(f = of, grad = og, hess = oh)
}

k <- 30
vprobs <- make_vprob(k)
bt <- rbeta(k + 1, 0.1, 0.1); bt <- bt/sum(bt)
(bt <- bt[-1])
ft <- vprobs %*% c(1-sum(bt), bt)
ys <- sample(0:k, 10000, TRUE, prob = ft)
(ws <- sapply(0:k, function(i) sum(ys == i)))
# fs <- make_functions(vprobs, ft, 0)
fs <- make_functions(vprobs, ws, 0.01)
fs$f(bt)
fs$f(rep(1/(k+1), k))

res_nlm <- suppressWarnings(nlm(fs$f, rep(1/(k+1), k)))
res_nlm$minimum
err1 <- l1d(bt, res_nlm$estimate)

res_nlm2 <- optim(rep(1/(k+1), k), fs$f)
res_nlm2
err2 <- l1d(bt, res_nlm2$par)

res_nlm3 <- optim(rep(1/(k+1), k), fs$f, gr = fs$grad, method = "BFGS")
res_nlm3
err3 <- l1d(bt, res_nlm3$par)

fs0 <- make_functions(vprobs, ws, 0)
res_nlm4 <- constrOptim(rep(1/(k+1), k), fs0$f, grad = fs0$grad,
                        ui = rbind(pracma::eye(k), -rep(1, k)),
                        ci = c(rep(0, k), -1), outer.eps = 1e-7)
res_nlm4$value
fs0$f(bt)
err4 <- l1d(bt, res_nlm4$par)


c(err1, err2, err3, err4)



####
##  Convert MLE back into original space
####


res_nlm <- suppressWarnings(nlm(fs$f, rep(1/(k+1), k)))
bth <- res_nlm$estimate
ft <- vprobs %*% c(1-sum(bth), bth)
bth <- c(1-sum(bth), bth)
dgu <- bth / ((k + 1):1)
gu <- cumsum(dgu)
sum(gu)
binprobs <- matrix(0, k + 1, k + 1)
for (i in 0:k) binprobs[, i + 1] <- dbinom(0:k, k, i/k)
cbind(ft, binprobs %*% gu)
