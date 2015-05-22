source("eb_ident/compare2.R")

## Test the code
run_exp(n= 30, p= 100, q= 2, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)
run_exp(n= 30, p= 200, q= 2, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)

run_exp(n= 100, p= 10, q= 10, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)
run_exp(n= 100, p= 20, q= 10, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)

library(parallel)
(res <- mclapply(1:7,
  function(i) {
    set.seed(i)
    run_exp(n= 50, p= 20, q= 20, W_X= 4, s_e= 2, W_e= 4, L= 100, n_te= 100)
  }, mc.cores = 7))