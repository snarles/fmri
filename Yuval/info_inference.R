####
##  Precision paper
##  identification loss with low-dim, high-dim, and mixture
####

library(pracma)
library(lineId)
source("info_theory_sims/fit_ident_curve.R")

## low-dim

load("Yuval/scores.RData")  ## obtained by running Yuval/identification2a.R


ks <- 20 * 1:12

infos_hd <- list()
infos_id <- list()

for (i in 1:5) {
  dmat <- scoreZ[[i]]
  nte <- 250
  #dmat <- pdist2(Yhat, Yte)
  acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
  (ihats_hd <- sapply(ks, function(k) Ihat_LI(1 - acs[k], k)))
  (ihats_id <- sapply(ks, function(k) aba_to_mi_lower(k, acs[k])))
  infos_hd[[i]] <- ihats_hd
  infos_id[[i]] <- ihats_id
  
  
  pdf(paste0("Yuval/ident_infer", i, ".pdf"))
  plot(ks, ihats_hd, ylim = c(0, 7), type = "o", xlim = c(0, 250), xlab = "k",
       ylab = "information")
  points(ks, ihats_id, type = "o", col = "blue", xlim = c(0, 250))
  title(paste(i * 100, "top voxels"))
  legend(170, 1, c("HD", "Ident"), col = c("black", "blue"), lwd = 1, pch = "o")
  dev.off()
  
  #plot(acs, type = "l")
  # k_sub <- 125
  # (I_implied <- fit_I_to_curve(acs[1:k_sub], wt_exp = 0))
  # acs_hat <- acs_curve(I_implied, 1:nte)
  # pdf(paste0("Yuval/ident_extrap", i, ".pdf"))
  # plot(1:nte, acs, type = "l", ylim = c(0, 1),
  #      xlab = "k", ylab = "accuracy")
  # legend(150, 0.5, col = c("black", "red"), lwd = 1, 
  #        legend = c("empirical", "extrapolated"))
  # lines(acs_hat, col = "red"); abline(v = k_sub)
  # m_err <- sqrt(mean(((acs_hat - acs)^2)[-(1:k_sub)]))
  # title(paste(i * 100, "top voxels"), 
  #       sub = paste("RMSE =", m_err))
  # dev.off()
}

infos_hd2 <- do.call(cbind, infos_hd)
infos_id2 <- do.call(cbind, infos_id)

pdf("Yuval/info_infer.pdf")
plot(100 * 1:5, infos_hd2[1, 1:5], ylim = c(0, 7), type = "o",
     xlab = "no. voxels", ylab = "information")
points(100 * 1:5, infos_id2[1, 1:5], col = "blue", type = "o")
legend(350, 1, c("HD", "Ident"), col = c("black", "blue"), lwd = 1, pch = "o")
dev.off()
