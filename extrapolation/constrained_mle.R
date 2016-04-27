make_vprob <- function(k, ps = seq(0, 1, 1/(2 * k))) {
  binprobs <- matrix(0, k + 1, length(ps))
  for (i in 1:length(ps)) binprobs[, i] <- dbinom(0:k, k, ps[i])
  # matplot(binprobs, type = "l")
  vprobs <- matrix(0, k + 1, length(ps))
  for (i in 1:length(ps))
    vprobs[, i] <- (binprobs %*% c(rep(0, i - 1),
                                   rep(1, length(ps) + 1 - i)))/(length(ps) + 1 - i)
  vprobs
}

cons_mle_est <- function(ppmat, k, ps = seq(0, 1, 1/(2 * k)), lbda = 0.1) {
  Ys <- as.numeric(ppmat)
  (ws <- sapply(0:k, function(i) sum(ys == i)))
  vprobs <- make_vprob(k, ps)
  of <- function(bt) {
    ft <- vprobs %*% c(1-sum(bt), bt)
    -sum(ws * log(ft)) - lbda * (sum(log(bt)) + log(1 - sum(bt)))
  }
  est <- suppressWarnings(nlm(of, rep(1/(length(ps)+1), length(ps) - 1)))
  bt <- est$estimate
  # ft <- vprobs %*% c(1-sum(bt), bt)
  bt <- c(1-sum(bt), bt)
  dgu <- bt / (length(ps):1)
  gu <- cumsum(dgu)
  # binprobs <- matrix(0, k + 1, length(ps))
  # for (i in 1:length(ps)) binprobs[, i] <- dbinom(0:k, k, ps[i])
  # ft2 <- binprobs %*% gu
  list(ps = ps, gu = gu)
}

cm_est_moment <- function(cm, K) {
  sum(cm$gu * cm$ps^K)
}

