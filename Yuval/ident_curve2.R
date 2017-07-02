####
##  Precision paper
##  identification loss with low-dim, high-dim, and mixture
####

library(pracma)
library(lineId)
source("info_theory_sims/fit_ident_curve.R")

## low-dim

load("Yuval/scores.RData")  ## obtained by running Yuval/identification2a.R

mugrid <- seq(2.5, 4, 0.05)
sigma2grid <- seq(0.8, 1.6, 0.05)
length(mugrid) * length(sigma2grid)


i <- 1
plot(1:250, acs, type = "l", ylim = c(0, 1))
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 2.9, sigma2 = 1), col = "red")
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 3.0, sigma2 = 1.5), col = "blue")
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 2.8, sigma2 = 0.7), col = "green")

i <- 2
plot(1:250, acs, type = "l", ylim = c(0, 1))
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 3.35, sigma2 = 1), col = "red")
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 3.45, sigma2 = 1.4), col = "blue")

i <- 3
plot(1:250, acs, type = "l", ylim = c(0, 1))
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 3.43, sigma2 = 1), col = "red")
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 3.6, sigma2 = 1.4), col = "blue")

i <- 4
plot(1:250, acs, type = "l", ylim = c(0, 1))
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 3.45, sigma2 = 1), col = "red")
lines(1:250, 1-piK(K = 1:250, mc.reps = 1000, mus = 3.65, sigma2 = 1.4), col = "blue")

for (i in 1:10) {
  dmat <- scoreZ[[i]]
  nte <- 250
  #dmat <- pdist2(Yhat, Yte)
  acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
  #plot(acs, type = "l")
  k_sub <- 125
  
  plot(1:250, acs, type = "l", ylim = c(0, 1))
  
  (I_implied <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
  acs_hat <- acs_curve(I_implied, 1:nte)
  pdf(paste0("Yuval/ident_extrap", i, ".pdf"))
  plot(1:nte, acs, type = "l", ylim = c(0, 1),
       xlab = "k", ylab = "accuracy")
  legend(150, 0.5, col = c("black", "red"), lwd = 1, 
         legend = c("empirical", "extrapolated"))
  lines(acs_hat, col = "red"); abline(v = k_sub)
  m_err <- sqrt(mean(((acs_hat - acs)^2)[-(1:k_sub)]))
  title(paste(i * 100, "top voxels"), 
        sub = paste("RMSE =", m_err))
  dev.off()
}
