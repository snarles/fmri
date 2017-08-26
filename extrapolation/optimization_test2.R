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
us <- seq(0, 1, 0.1)
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

library(nloptr)
gof <- function(gu) {
  ft <- binprobs %*% gu
  -as.numeric(t(binprobs) %*% (ws/ft))
}
of_gu(gu_unif)
gof(gu_unif)
numDeriv::grad(of_gu, gu_unif)
res <- nloptr(gu_unif, of_gu, eval_grad_f = gof, 
              lb = 0 * us, ub = 0 * us + 1, 
              eval_g_ineq = function(gu) sum(gu) - 1, 
              eval_jac_g_ineq = function(gu) 0 * gu + 1, 
              opts = list(algorithm = "NLOPT_LD_MMA",
                          xtol_rel = 1.0e-8,
                          print_level = 2,
                          check_derivatives = TRUE,
                          check_derivatives_print = "all"))
print(res)
gu_temp <- res$solution
of_gu(gu_temp)
of_gu(gu)
plot(gu_temp)
lines(gu)


# library(alabama)
# hinf <- function(gu) gu
# heqf <- function(gu) 1 - sum(gu)
# gu_temp <- gu_unif
# res <- auglag(gu_temp, of_gu, hin = hinf, heq = heqf)
# 
# library(NlcOptim)
# library(pracma)
# gu_temp <- gu_unif2
# res <-NlcOptim(gu_temp, of_gu, confun = function(x) return(list(ceq=0, c=1)),
#                A = -eye(length(us)),
#                B = t(t(0 * us)), Aeq = t(0 * us + 1), Beq = 1)


# library(Rsolnp)
# gu_temp <- gu_unif
# res <- solnp(gu_temp, of_gu, eqfun = heqf, eqB = 1, LB = 0 * us, UB = 1 + 0 * us)


# library(nloptr)
# gu_temp <- gu_unif
# res <- auglag(gu_temp, of_gu, hin = hinf, heq = heqf)


####
###  TEST NLOPTR
####

library(nloptr)
eval_f <- function(x) {
  return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}
eval_grad_f <- function(x) {
  c(-400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
    200 * (x[2] - x[1]^2))
}
x0 <- c(-1.2, 1)
opts = list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1.0e-8)
res <- nloptr(x0 = x0, eval_f = eval_f, eval_grad_f = eval_grad_f, opts = opts)

## Second example

eval_f0 <- function(x, a, b) {
  return( sqrt(x[2]) )
}
eval_grad_f0 <- function(x, a, b) {
  return( c(0, .5/sqrt(x[2])) )
}
eval_g0 <- function(x, a, b) {
  return( (a * x[1] + b)^3 - x[2] )
}
eval_jac_g0 <- function(x, a, b) {
  return( rbind( c(3*a[1]*(a[1]*x[1] + b[1])^2, -1.0 ),
                 c(3*a[2]*(a[2]*x[2] + b[2])^2, -1.0 ) ) )
}
a <- c(2, -1)
b <- c(0, 1)
res0 <- nloptr(x0 = c(1.234, 5.678),
               eval_f = eval_f0,
               eval_grad_f = eval_grad_f0,
               lb = c(-Inf, 0),
               ub = c(Inf, Inf),
               eval_g_ineq = eval_g0,
               eval_jac_g_ineq = eval_jac_g0, 
               opts = list(algorithm = "NLOPT_LD_MMA",
                           xtol_rel = 1.0e-8,
                           print_level = 2,
                           check_derivatives = TRUE,
                           check_derivatives_print = "all"),
               a = a, b = b)


