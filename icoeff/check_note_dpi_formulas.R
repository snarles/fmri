## check that 
## d/de D_f(Q^e p||Q^e q) = \int [f(p/q) - (p/q)f'(p/q)] (Qq-q) + f'(p/q) (Qp - p)  dx

f <- function(x) {
  ans <- x * log(x)
  ans[x==0] <- 0
  ans
}
f1 <- function(x) 1 + log(x)

sim_delt <- 1e-3
xs <- seq(-5, 5, sim_delt)

f_div <- function(dens_p, dens_q) {
  sim_delt * sum(f(dens_p/dens_q) * dens_q)  
}

formula_rhs <- function(dens_p, dens_q, dens_Qp, dens_Qq) {
  sim_delt * sum(
    (f(dens_p/dens_q) - dens_p/dens_q * f1(dens_p/dens_q)) * (dens_Qq - dens_q) + 
      f1(dens_p/dens_q) * (dens_Qp - dens_p)
  )
}

formula_rhs2 <- function(dens_p, dens_q, dens_Qp, dens_Qq) {
  sim_delt * sum(
    f(dens_p/dens_q) * (dens_Qq - dens_q) + 
      f1(dens_p/dens_q) * (dens_Qp - (dens_p/dens_q)*dens_Qq)
  )
}

formula_rhs_sp <- function(dens_p, dens_q, dens_Qp, dens_Qq) {
  sim_delt * sum(
    (dens_Qp - dens_p) * log(dens_p/dens_q) + dens_Qp - (dens_p/dens_q)*dens_Qq
  )
}


mu1 <- 0; sigma2_1 <- 1
mu2 <- 0.2; sigma2_2 <- 1

sigma2_Q <- 0.1

dens_p <- dnorm(xs, mean = mu1, sd = sqrt(sigma2_1))
dens_q <- dnorm(xs, mean = mu2, sd = sqrt(sigma2_2))
dens_Qp <- dnorm(xs, mean = mu1, sd = sqrt(sigma2_1 + sigma2_Q))
dens_Qq <- dnorm(xs, mean = mu2, sd = sqrt(sigma2_2 + sigma2_Q))

num_deltas <- 10^(-(1:10))

res <- sapply(num_deltas, function(num_delta) {
  dens_Qdp <- (1-num_delta) * dens_p + num_delta * dens_Qp
  dens_Qdq <- (1-num_delta) * dens_q + num_delta * dens_Qq
  
  f_div(dens_p, dens_q)
  f_div(dens_Qdp, dens_Qdq)
  c((f_div(dens_Qdp, dens_Qdq) - f_div(dens_p, dens_q))/num_delta,
    formula_rhs(dens_p, dens_q, dens_Qp, dens_Qq),
    formula_rhs2(dens_p, dens_q, dens_Qp, dens_Qq),
    formula_rhs_sp(dens_p, dens_q, dens_Qp, dens_Qq))
})

matplot(-log(num_deltas), t(res), type = "l")
