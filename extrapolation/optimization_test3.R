library(nloptr)

binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

## Generate data
n <- 1000
us <- seq(0, 1, 0.1)
gu <- cumsum(runif(length(us)) * (1:length(us))^4)
gu <- gu/sum(gu)
plot(us, gu, type = "l")
usamp <- sample(us, n, TRUE, prob = gu)
k <- 20
Ys <- rbinom(n, k, usamp)
(ws <- sapply(0:k, function(i) sum(Ys == i)))
(momK <- mean(binmom(Ys, k, k)))
usk <- us^k

disp_solution <- function(gu_temp) {
  print(c(new = of_gu(gu_temp), true = of_gu(gu)))
  plot(gu_temp)
  lines(gu)
}

## Setup objective function
binprobs <- matrix(0, k + 1, length(us))
for (i in 1:length(us)) binprobs[, i] <- dbinom(0:k, k, us[i])
of_gu <- function(gu) {
  ft <- binprobs %*% gu
  -sum(ws * log(ft))
}
gof <- function(gu) {
  ft <- binprobs %*% gu
  -as.numeric(t(binprobs) %*% (ws/ft))
}

## Initial objective value
gu_unif <- rep(1/length(us), length(us))
of_gu(gu)
of_gu(gu_unif)


## MPLE unconstrained
res <- nloptr(gu_unif, of_gu, eval_grad_f = gof, 
              lb = 0 * us, ub = 0 * us + 1, 
              eval_g_ineq = function(gu) sum(gu) - 1, 
              eval_jac_g_ineq = function(gu) 0 * gu + 1, 
              opts = list(algorithm = "NLOPT_LD_MMA",
                          xtol_rel = 1.0e-8,
                          print_level = 2,
                          check_derivatives = TRUE,
                          check_derivatives_print = "all"))
# print(res)
gu_temp <- res$solution
disp_solution(res$solution)


## MPLE moment constraint
res <- nloptr(gu_unif, of_gu, eval_grad_f = gof, 
              lb = 0 * us, ub = 0 * us + 1, 
              eval_g_ineq = function(gu) sum(gu) - 1, 
              eval_jac_g_ineq = function(gu) 0 * gu + 1, 
              eval_g_eq = function(gu) sum(gu * usk) - momK,
              eval_jac_g_eq = function(gu) usk,
              opts = list(algorithm = "NLOPT_LD_MMA",
                          xtol_rel = 1.0e-8,
                          print_level = 2,
                          check_derivatives = TRUE,
                          check_derivatives_print = "all"))
# print(res)
gu_temp <- res$solution
disp_solution(res$solution)