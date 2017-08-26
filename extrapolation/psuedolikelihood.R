library(nloptr)

binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}

fit_pm_models <- function(Ys, k, us = seq(0, 1, 0.02), gu_init = rep(1/length(us), length(us))) {
  Ys <- as.numeric(Ys)
  (ws <- sapply(0:k, function(i) sum(Ys == i)))
  (momK <- mean(binmom(Ys, k, k)))
  usk <- us^k
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
  ## MPLE unconstrained
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = function(gu) sum(gu) - 1, 
                eval_jac_g_ineq = function(gu) 0 * gu + 1, 
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,
                            check_derivatives = FALSE,maxeval = 1e4))
  (t2u <- proc.time() - t1)
  # print(res)
  gu_mple <- res$solution
  # disp_solution(res$solution); title("uncon")
  
  
  ## MPLE moment constraint
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = function(gu) sum(gu) - 1, 
                eval_jac_g_ineq = function(gu) 0 * gu + 1, 
                eval_g_eq = function(gu) sum(gu * usk) - momK,
                eval_jac_g_eq = function(gu) usk,
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,
                            check_derivatives = FALSE,maxeval = 1e4))
  (t2mom <- proc.time() - t1)
  # print(res)
  gu_mom <- res$solution
  
  ## MPLE monotonic constraint
  ll <- length(us)
  mat <- matrix(0, ll - 1, ll)
  cmat <- (row(mat) == col(mat)) - (row(mat) == (col(mat) - 1))
  eval_g_ineq_mono = function(gu) c(sum(gu) - 1, gu[-ll] - gu[-1])
  eval_jac_g_ineq_mono = function(gu) rbind(0 * gu + 1, cmat)
  
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = eval_g_ineq_mono, 
                eval_jac_g_ineq = eval_jac_g_ineq_mono, 
                # eval_g_eq = function(gu) sum(gu * usk) - momK,
                # eval_jac_g_eq = function(gu) usk,
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,
                            check_derivatives = FALSE,maxeval = 1e4))
  (t2mono <- proc.time() - t1)
  # print(res)
  gu_mono <- res$solution

  ## MPLE 2 constraint
  ll <- length(us)
  mat <- matrix(0, ll - 1, ll)
  cmat <- (row(mat) == col(mat)) - (row(mat) == (col(mat) - 1))
  eval_g_ineq_mono = function(gu) c(sum(gu) - 1, gu[-ll] - gu[-1])
  eval_jac_g_ineq_mono = function(gu) rbind(0 * gu + 1, cmat)
  
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = eval_g_ineq_mono, 
                eval_jac_g_ineq = eval_jac_g_ineq_mono, 
                eval_g_eq = function(gu) sum(gu * usk) - momK,
                eval_jac_g_eq = function(gu) usk,
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,check_derivatives = FALSE,
                            maxeval = 1e4))
  (t2c2 <- proc.time() - t1)
  # print(res)
  gu_mm <- res$solution
  
  list(gu_mple = gu_mple, gu_mom = gu_mom, gu_mono = gu_mono, gu_mm = gu_mm)
}