source("eb_ident/compare2.R")

#pars <- list(n= 30, p= 100, q= 5, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)
#attach(pars)
#detach(pars)

## Test the code

run_exp(n= 30, p= 100, q= 2, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)
run_exp(n= 30, p= 200, q= 2, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)

run_exp(n= 100, p= 10, q= 10, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)
run_exp(n= 100, p= 20, q= 10, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)

# large-scale experiments

library(parallel)

def_settings <- function(xx) {
  set.seed(xx)
  run_exp(n=driver_n, p= 60, q= 60, W_X= 2, s_e= 10, W_e= 2, L= 100, n_te= 100)  
}

n_seeds <- 7
n_iters <- 10
results <- array(0, dim = c(n_seeds, 5, n_iters))

for (ii in 1:n_iters) {
  driver_n <- ii * 10
  res <- mclapply(1:n_seeds, def_settings, mc.cores = 7)
  results[, , ii] <- t(matrix(unlist(res), nrow = 5))  
}

