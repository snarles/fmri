####
##  Precision paper
##  identification loss with low-dim, high-dim, and mixture
####

library(pracma)
library(lineId)
library(parallel)
source("info_theory_sims/fit_ident_curve.R")

## low-dim

load("Yuval/scores_MSE_c.RData")  ## obtained by running Yuval/identification2a.R

mugrid <- seq(2.7, 6, 0.01)
sigma2grid <- seq(1, 5, 0.025)
length(mugrid) * length(sigma2grid)

parmat <- cbind(mu = rep(mugrid, each = length(sigma2grid)),
                sigma2 = rep(sigma2grid, length(mugrid)))

nrow(parmat)

t1 <- proc.time()
par_res <- mclapply(1:nrow(parmat), function(i) {
  mu <- parmat[i, 1]
  sigma2 <- parmat[i, 2]
  1-piK(K = 1:250, mc.reps = 1000, mu, sigma2)
}, mc.cores = 2)
proc.time() - t1

precomputed_curves <- do.call(cbind, par_res)
save(parmat, precomputed_curves, file = "Yuval/ident_curve_precomp.rda")




fitted_pars <- matrix(0, 10, 2)
fitted_mus <- numeric(10)

for (i in 1:10) {
  dmat <- scoreZ[[i]]
  nte <- 250
  #dmat <- pdist2(Yhat, Yte)
  acs <- 1-sapply(1:nte, function(i) resample_misclassification(-dmat, 1:nte, i))
  #plot(acs, type = "l")
  k_sub <- 125
  
  ## find best 2-parameter fit
  dd <- colSums((precomputed_curves - acs)[1:k_sub, ]^2)
  pars <- parmat[which.min(dd), ]
  acs_hat_2p <- precomputed_curves[, which.min(dd)]
  
  ## find best fit with sigma=1
  dd[parmat[, "sigma2"] != 1] <- Inf
  mu_fit <- parmat[which.min(dd), 1]
  acs_hat_1p <- precomputed_curves[, which.min(dd)]

  #pdf(paste0("Yuval/ident_extrap_ip", i, ".pdf"))
  plot(1:250, acs, type = "l", ylim = c(0, 1))
  lines(acs_hat_2p, col = "blue")
  lines(acs_hat_1p, col = "red")
  abline(v = k_sub)
  legend(150, 0.5, col = c("black", "red", "blue"), lwd = 1, 
         legend = c("empirical", "1par", "2par"))
  title(paste(i * 100, "top voxels"))
  #dev.off()
  
  fitted_pars[i, ] <- pars
  fitted_mus[i] <- mu_fit
}

cbind(mu_1 = fitted_mus, mu_2 = fitted_pars[, 1], sigma2 = fitted_pars[, 2])
