source("extrapolation/fit_expo.R")

## binomial moment fomula
binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

mle_est <- function(ppmat, k, ps = seq(0, 1, 1/(2 * k))) {
  Ys <- as.numeric(ppmat)
  (ws <- sapply(0:k, function(i) sum(Ys == i)))
  binprobs <- matrix(0, k + 1, length(ps))
  for (i in 1:length(ps)) binprobs[, i] <- dbinom(0:k, k, ps[i])
  of_gu <- function(gu) {
    ft <- binprobs %*% gu
    -sum(ws * log(ft))
  }
  
  bt <- est$estimate
  # ft <- vprobs %*% c(1-sum(bt), bt)
  bt <- c(1-sum(bt), bt)
  dgu <- bt / (length(ps):1)
  gu <- cumsum(dgu)
  # ft2 <- binprobs %*% gu
  list(ps = ps, gu = gu, of = of, of_gu = of_gu)
}

n <- 1000
us <- seq(0, 1, 0.05)
gu <- cumsum(runif(length(us)) * (1:length(us))^4)
gu <- gu/sum(gu)
plot(us, gu, type = "l")
usamp <- sample(us, n, TRUE, prob = gu)
k <- 20
Ys <- rbinom(n, k, usamp)
(ws <- sapply(0:k, function(i) sum(Ys == i)))

binprobs <- matrix(0, k + 1, length(us))
for (i in 1:length(us)) binprobs[, i] <- dbinom(0:k, k, us[i])
of_gu <- function(gu) {
  ft <- binprobs %*% gu
  -sum(ws * log(ft))
}
penalty_pos <- 0.1; penalty_one <- 0.1

scl <- 1
of_gu_pen <- function(gu) {
  ft <- binprobs %*% gu
  -scl * sum(ws * log(ft)) + penalty_pos * sum(-log(gu)) + penalty_one * sum(-log(1 - sum(gu)))
}

gu_unif <- rep(1/length(us), length(us))
of_gu(gu)
of_gu(gu_unif)
gu_unif2 <- rep(1/length(us), length(us)) * 0.999
of_gu_pen(gu_unif2)

help(optim)
gu_temp <- gu_unif2
gu_temp <- optim(gu_temp, of_gu_pen, method = "Nelder-Mead", control = list(trace = 1,
                                                                     maxit = 1e5))$par

of_gu(gu_temp)
of_gu(gu)

penalty_one <- 1e-4; penalty_pos <- 1e-4
penalty_one <- 1e-5; penalty_pos <- 1e-5
penalty_one <- 1e-6; penalty_pos <- 1e-6
scl <- 3e3

of_gu_pen(gu * 0.999)
plot(gu_temp)
penalty_pos <- 0.01; penalty_one <- 0.01
gu_temp <- optim(gu_temp, of_gu_pen)$par
of_gu(gu_temp)
plot(gu_temp)
penalty_pos <- 1e-4; penalty_one <- 1e-4
gu_temp <- optim(gu_temp, of_gu_pen, method = "SANN", control = list(trace = 1))$par
of_gu(gu_temp)
of_gu_pen(gu_temp)

plot(gu_temp)
lines(gu)
