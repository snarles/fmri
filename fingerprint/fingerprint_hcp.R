library(pracma)
library(MASS)
library(class)
library(parallel)
library(lineId) ## use devtools::install('lineId')
source('extrapolation/kay_method2.R')


gsr_kl <- -as.matrix(read.csv('fingerprint/FL_GSR_kl.csv'))
gsr_cor <- as.matrix(read.csv('fingerprint/FL_GSR_cor.csv'))
ngsr_kl <- -as.matrix(read.csv('fingerprint/FL_noGSR_kl.csv'))
ngsr_cor <- as.matrix(read.csv('fingerprint/FL_noGSR_cor.csv'))


nsubs <- dim(gsr_kl)[1] # 338




mi_est_pipeline <- function(pmat) {
  K <- dim(pmat)[1]
  empirical_acc <- 1- resample_misclassification(pmat, 1:K)
  
  
  acc_crit <- 0.75
  # extrapolate until 0.75 acc
  kgrid <- floor(exp((1:300)/5))
  accs_kay <- kernel_extrap3(pmat, kgrid)
  #accs_kay
  #plot(accs_kay)
  
  kchosen <- min(kgrid[accs_kay < acc_crit])
  acc_at_k <- max(accs_kay[accs_kay < acc_crit])
  
  ## get MI estimate
  #mi_est <- Ihat_LI(1-acc_at_k, kchosen)
  mi_grid <- log(kchosen) * (1:100)/200
  pi_ests <- piK(mi_grid, kchosen)
  mi_est <- max(mi_grid[(1-pi_ests) < acc_crit ])
  
  list(empirical_acc=empirical_acc, kchosen=kchosen, acc_at_k=acc_at_k, 
       mi_est = mi_est, acc_at_est = 1-piK(mi_est, kchosen))
}


est_gsr_kl <- numeric()
est_ngsr_kl <- numeric()
est_gsr_cor <- numeric()
est_ngsr_cor <- numeric()

nreps <- 100
nsubsample <- 20
all_inds <- matrix(nrow=nreps, ncol = nsubsample)

variable <- 'empirical_acc'
#variable <- 'mi_est'


for (i in 1:nreps) {
  # scores matrix
  set.seed(i)
  inds <- sample(nsubs, nsubsample)
  all_inds[i, ] <- inds
  est_gsr_kl[i] <- mi_est_pipeline(gsr_kl[inds, inds])[[variable]]
  est_ngsr_kl[i] <- mi_est_pipeline(ngsr_kl[inds, inds])[[variable]]
  est_gsr_cor[i] <- mi_est_pipeline(gsr_cor[inds, inds])[[variable]]
  est_ngsr_cor[i] <- mi_est_pipeline(ngsr_cor[inds, inds])[[variable]]
}

results <- cbind(est_gsr_kl, est_ngsr_kl, est_gsr_cor, est_ngsr_cor)
hist(results[,1] - results[,2])
hist(results[,3] - results[,4])

