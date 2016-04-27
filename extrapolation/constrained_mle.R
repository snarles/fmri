make_vprob <- function(k) {
  binprobs <- matrix(0, k + 1, k + 1)
  for (i in 0:k) binprobs[, i + 1] <- dbinom(0:k, k, i/k)
  # matplot(binprobs, type = "l")
  vprobs <- matrix(0, k + 1, k + 1)
  for (i in 0:k) vprobs[, i+1] <- (binprobs %*% c(rep(0, i), rep(1, k + 1 - i)))/(k + 1 - i)
  vprobs
}

cons_mle_est <- function(ppmat, k, lambda) {
  Ys <- as.numeric(ppmat)
  (ws <- sapply(0:k, function(i) sum(ys == i)))
  vprobs <- make_vprob(k)
  of <- function(bt) {
    ft <- vprobs %*% c(1-sum(bt), bt)
    -sum(ws * log(ft)) - lambda * (sum(log(bt)) + log(1 - sum(bt)))
  }
  est <- suppressWarnings(nlm(of, rep(1/(k+1), k)))
  bt <- est$estimate
  bt <- c(1-sum(bt), bt)
  dgu <- bt / ((k + 1):1)
  gu <- cumsum(dgu)
  gu
}

cm_est_moment <- function(gu, K) {
  k <- length(gu) - 1
  ps <- (0:k)/k
  sum(gu * ps^K)
}

