library(pracma)
library(lineId)
source('extrapolation/kay_method2.R')


mi_est_pipeline <- function(pmat) {
  K <- dim(pmat)[1]
  empirical_acc <- 1- resample_misclassification(pmat, 1:K)
  
  
  acc_crit <- 0.75
  # extrapolate until 0.75 acc
  kgrid <- floor(exp((1:600)/20))
  accs_kay <- kernel_extrap3(pmat, kgrid)
  #accs_kay
  #plot(accs_kay)
  
  kchosen <- min(kgrid[accs_kay < acc_crit])
  acc_at_k <- max(accs_kay[accs_kay < acc_crit])
  
  ## get MI estimate
  #mi_est <- Ihat_LI(1-acc_at_k, kchosen)
  mi_grid <- log(kchosen) * (1:300)/200
  pi_ests <- piK(sqrt(2 * mi_grid), kchosen)
  mi_est <- max(mi_grid[(1-pi_ests) < acc_crit ])
  
  list(empirical_acc=empirical_acc, kchosen=kchosen, acc_at_k=acc_at_k, 
       mi_est = mi_est, acc_at_est = 1-piK(sqrt(2 * mi_est), kchosen))
}


mi_est_pipeline2 <- function(pmat) {
  K <- dim(pmat)[1]
  empirical_acc <- 1- resample_misclassification(pmat, 1:K)
  accs_kay <- kernel_extrap3(pmat, K)

  ## get MI estimate
  #mi_est <- Ihat_LI(1-acc_at_k, kchosen)
  mi_grid <- log(K) * (1:800)/200
  pi_ests <- piK(sqrt(2 * mi_grid), K)
  mi_est <- max(mi_grid[(1-pi_ests) < accs_kay ])
  
  list(empirical_acc=empirical_acc, accs_kay=accs_kay, 
       mi_est = mi_est, acc_at_est = 1-piK(sqrt(2 * mi_est), K), max_mi = max(mi_grid))
}