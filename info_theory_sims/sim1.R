####
##  ESTIMATING MI FROM BAYES ERROR/BAYES CONFUSION
####

source("info_theory_sims/gaussian_source.R")
source("info_theory_sims/methods_source.R")

####
##  Gaussian
####

d <- 10
sigma2 <- 5
K <- 5

(i_true <- d/2 * log(1 + (1/sigma2)))
(snr <- d/2/sigma2)
(abe <- mc_ident2(d, sigma2, K))
(ihat_LS <- Ihat_LS(abe, K))
fK(sqrt(2 * ihat_LS), K)
(ihat_fano <- Ihat_fano(abe, K))



