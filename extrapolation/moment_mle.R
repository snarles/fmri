####
##  Moment-constrained MLE
#### 

binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

momk_mle_est <- function(ppmat, k, ps = seq(0, 1, 1/(2 * k)), lbda = 0.1, mpen = 10) {
  Ys <- as.numeric(ppmat)
  momk <- mean(binmom(Ys, k, k))
  psk <- ps^k
  (ws <- sapply(0:k, function(i) sum(Ys == i)))
  binprobs <- matrix(0, k + 1, length(ps))
  for (i in 1:length(ps)) binprobs[, i] <- dbinom(0:k, k, ps[i])
  of <- function(gum) {
    gu <- c(gum, 1-sum(gum))
    guk <- sum(psk * gu)
    ft <- binprobs %*% gu
    -sum(ws * log(ft)) - lbda * sum(log(gu)) + 
      mpen * (guk - momk)^2
  }
  gum <- rep(1/(length(ps)), length(ps) - 1)
  est <- suppressWarnings(nlm(of, gum))
  gum <- est$estimate
  gu <- c(gum, 1-sum(gum))
  sum(gu * psk)
  list(gu = gu, ps = ps, momk = c(momk, sum(gu * psk)))
}

