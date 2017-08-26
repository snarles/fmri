####
##  Precision paper
##  identification loss with low-dim, high-dim, and mixture
####

library(pracma)
library(lineId)
source("info_theory_sims/fit_ident_curve.R")

## low-dim

load("Yuval/scores.RData")  ## obtained by running Yuval/identification2a.R

for (i in 1:10) {
  dmat <- scoreZ[[i]]
  nte <- 250
  #dmat <- pdist2(Yhat, Yte)
  acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
  #plot(acs, type = "l")
  k_sub <- 125
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
