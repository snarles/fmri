####
##  ESTIMATING MI FROM BAYES ERROR/BAYES CONFUSION
####

source("info_theory_sims/gaussian_source.R")
source("info_theory_sims/methods_source.R")
library(parallel); mcc <- 3


####
##  Gaussian
####



d <- 10
sigma2 <- 0.5
K <- 5
cm.reps <- 100

(i_true <- d/2 * log(1 + (1/sigma2)))
(snr <- d/2/sigma2)
(abe <- mc_ident2(d, sigma2, K))
(ihat_LS <- Ihat_LS(abe, K))
fK(sqrt(2 * ihat_LS), K)
(ihat_fano <- Ihat_fano(abe, K))
ihat_cm_temp <- unlist(mclapply(1:cm.reps, 
                         function(i) {
                           cm <- mc_cm(d, sigma2, K)
                           cm_to_I(cm)
                         }, mc.cores = mcc))
(ihat_cm <- mean(ihat_cm_temp))

c(i_true = i_true, ihat_LS = ihat_LS, ihat_fano = ihat_fano, ihat_cm = ihat_cm)

