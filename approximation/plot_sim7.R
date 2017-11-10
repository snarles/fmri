
source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")
source("par2/objective_function.R")
source("extrapolation/basis_source.R")

p <- 10
sigma2_seq <- 0.01 * 1:50
K <- 10000 ## multiple of 1000
Ktarg <- c(1000, 2000, 5000, 10000)
ksub <- 500 ## multiple of 250
ksub2 <- ksub/2

nboot <- 20 ## number of bootstraps for CV
mc.reps <- 100000
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))

kseq <- function(nr, ksub) {
  interv <- floor(ksub/nr)
  seq(interv, ksub, by = interv)
}

nspline <- 100
nrow <- 100
kref <- kseq(nrow, ksub)
kde_bdwids <- list("bcv", "ucv", 0.1, 0.15, 0.2, 0.25, 0.3)
(max.mu <- (qnorm(1- 1/(max(kref)^2))))
## fixed bandwidths
fixed.gb <- c(0.4, 0.45, 0.5, 0.55, 0.6)
lsub2 <- length(kref)/2
basis_vecs <- list()

## spline100
# {
#   ns <- nspline
#   nm <- paste0("r.lin",ns)
#   knts <- seq(0, 1, length.out = ns + 2)
#   knts <- rev(1 - knts[-c(1, ns + 2)]^2)
#   Xmat <- spline1_moments(knts, kref)
#   Xtarg <- spline1_moments(knts, Ktarg)
#   basis_vecs[[nm]] <- list(Xmat = Xmat, Xtarg = Xtarg)
# }

for (gb in fixed.gb) {
  nm <- paste0("r.gauss", gb)
  basis_vecs[[nm]] <- 0
}
for (gb in fixed.gb) {
  nm <- paste0("r.gauss", gb)
  basis_vecs[[nm]] <- get_basis_mat(max.mu, kernel_sigma = gb)
}

(bdwids <- seq(0.1, 1, 0.1))
# basis_sets <- lapply(bdwids, function(bd) get_basis_mat(max.mu, bd))
# 
# sub_basis_sets <- lapply(basis_sets, function(set1) {
#   list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(kref), , drop = FALSE])
# })

(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(names(basis_vecs), kde.names, "r.cv.gauss", "par2"))


res <- readRDS("approximation/raw_7_results.rds")


sigma2s <- sigma2s[1:length(res)]

resa <- do.call(rbind, lapply(res, `[[`, "accs"))
true_accs <- matrix(NA, length(sigma2_seq), ncol(resa))
for (ii in 1:length(sigma2_seq)) {
  true_accs[ii, ] <- colMeans(resa[sigma2s == sigma2_seq[ii], ])
}
matplot(Ktarg, t(true_accs), type = "l", ylim = c(0, 1))


## bootstrap the rmseZ


(facc <- true_accs)

rmseZ <- list()
biaseZ <- list()
all_final_predZ <- list()
rmse_sdZ <- list()
bias_sdZ <- list()
nboot <- 1000

for (ind in 1:length(Ktarg)) {
  all_final_preds <- do.call(rbind, lapply(res, function(v) v$preds[, ind]))
  #all_final_preds <- do.call(rbind, lapply(res, `[[`, "preds"))
  sum(all_final_preds < 0)
  all_final_preds[all_final_preds < 0] <- 0
  colnames(all_final_preds) <- column_names
  all_final_predZ[[ind]] <- all_final_preds
  resids <- all_final_preds - facc[,ind][match(sigma2s, sigma2_seq)]
  rmses <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
  biases <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
  rmse_sd <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
  bias_sd <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
  for (ii in 1:length(sigma2_seq)) {
    rmses[ii, ] <- sqrt(colMeans(resids[sigma2s == sigma2_seq[ii], ]^2))
    biases[ii, ] <- colMeans(resids[sigma2s == sigma2_seq[ii], ])
    rmse_boots <- matrix(NA, nboot, ncol(rmses))
    bias_boots <- matrix(NA, nboot, ncol(rmses))
    orig_inds <- which(sigma2s == sigma2_seq[ii])
    for (jj in 1:nboot) {
      boot_inds <- sample(orig_inds, length(orig_inds), replace = TRUE)
      rmse_boots[jj, ] <- sqrt(colMeans(resids[boot_inds, ]^2))
      bias_boots[jj,] <- colMeans(resids[boot_inds, ])
    }
    rmse_sd[ii, ] <- apply(rmse_boots, 2, sd)
    bias_sd[ii,] <- apply(bias_boots,2,sd)
  }
  colnames(rmses) <- column_names
  colnames(rmse_sd) <- column_names
  colnames(bias_sd) <- column_names
  colnames(biases) <- column_names
  rmseZ[[ind]] <- rmses
  rmse_sdZ[[ind]] <- rmse_sd
  biaseZ[[ind]] <- biases
  bias_sdZ[[ind]] <- bias_sd
}


minimax_rmses <- list()
minimax_rmse_sds <- list()
master_inds <- rep(1:200, each = 50)
nboot <- 1000

set.seed(0)
for (ind in 1:length(Ktarg)) {
  rmses <- rmseZ[[ind]]
  minimax_rmses[[ind]] <- apply(rmses, 2, max)
  resids <- all_final_predZ[[ind]] - facc[,ind][match(sigma2s, sigma2_seq)]
  mm_boot <- matrix(NA, nboot, ncol(all_final_preds))
  colnames(mm_boot) <- column_names
  for (ind.boot in 1:nboot) {
    boot.inds <- sample(200, 200, TRUE)
    rmses <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
    for (ii in 1:length(sigma2_seq)) {
      rmses[ii, ] <- sqrt(colMeans(resids[sigma2s == sigma2_seq[ii], ][boot.inds, ]^2))
    }
    mm_boot[ind.boot, ] <- apply(rmses, 2, max)
  }
  minimax_rmse_sds[[ind]] <- apply(mm_boot, 2, sd)
}

ind <- 1
ind <- 2
ind <- 3
ind <- 4
sel_vars <- c("r.cv.gauss", "kde_bcv", "kde_ucv")

rbind(minimax_rmses[[ind]][sel_vars], minimax_rmse_sds[[ind]][sel_vars])
# r.cv.gauss     kde_bcv     kde_ucv
# [1,] 0.032012576 0.090849604 0.067092107
# [2,] 0.001004484 0.001093817 0.001147322
# r.cv.gauss     kde_bcv     kde_ucv
# [1,] 0.043983286 0.087769100 0.059097824
# [2,] 0.001591107 0.001130881 0.001029907
# r.cv.gauss     kde_bcv      kde_ucv
# [1,] 0.073465930 0.079319514 0.0510825724
# [2,] 0.003623713 0.001112792 0.0009787097
# r.cv.gauss     kde_bcv      kde_ucv
# [1,] 0.097685388 0.075615097 0.0454476003
# [2,] 0.004153768 0.001218799 0.0009587312

library(reshape2)
library(ggplot2)

ind <- 1
ind <- 2
ind <- 3
ind <- 4
sd_mult <- 2.95
rmses <- rmseZ[[ind]]
sds <- rmse_sdZ[[ind]]
sel_vars <- c("r.cv.gauss", "kde_bcv", "kde_ucv")
temp <- data.frame(true_acc = true_accs[, ind], rmses[, sel_vars])
temp2 <- melt(data = temp, id.vars = "true_acc")
temp_se <- data.frame(true_acc = true_accs[, ind], sds[, sel_vars])
temp_se <- melt(data = temp_se, id.vars = "true_acc")
temp3 <- data.frame(temp2, rmse_low = temp2[, 3] - temp_se[, 3], rmse_high = temp2[, 3] + temp_se[, 3])
colnames(temp3)[3] <- "rmse"

ggplot(data = temp3, aes(x = true_acc, y = rmse, colour = variable, linetype=variable)) +
  geom_line() + coord_cartesian(xlim = c(0, 1)) + 
  geom_errorbar(aes(ymin = rmse_low, max = rmse_high)) + 
  scale_linetype_manual(values = c(1,4,6))+
  ggtitle(paste0("Predicting K=", Ktarg[ind], " from k=", ksub)) + ylim(c(0,0.11))
  #ggsave("approximation/sim_large7_K1_k0_5.png", width = 6, height = 3)
  #ggsave("approximation/sim_large7_K2_k0_5.png", width = 6, height = 3)
  #ggsave("approximation/sim_large7_K5_k0_5.png", width = 6, height = 3)
  #ggsave("approximation/sim_large7_K10_k0_5.png", width = 6, height = 3)
  

## plot biases type 2
ind <- 1
ind <- 2
ind <- 3
ind <- 4
sd_mult <- 2.95
biases <- biaseZ[[ind]]
sds <- bias_sdZ[[ind]]
sel_vars <- c("r.cv.gauss", "kde_bcv", "kde_ucv")
temp <- data.frame(true_acc = true_accs[, ind], biases[, sel_vars])
temp2 <- melt(data = temp, id.vars = "true_acc")
temp_se <- data.frame(true_acc = true_accs[, ind], sds[, sel_vars])
temp_se <- melt(data = temp_se, id.vars = "true_acc")
temp3 <- data.frame(temp2, bias_low = temp2[, 3] - temp_se[, 3], bias_high = temp2[, 3] + temp_se[, 3])
colnames(temp3)[3] <- "pred_minus_mean"

ggplot(data = temp3, aes(x = true_acc, y = pred_minus_mean, colour = variable, linetype=variable)) +
  geom_line() + coord_cartesian(xlim = c(0, 1)) + 
  geom_errorbar(aes(ymin = bias_low, max = bias_high)) + 
  scale_linetype_manual(values = c(1,4,6))+
  ggtitle(paste0("Predicting K=", Ktarg[ind], " from k=", ksub)) + ylim(c(-0.10,0.10))

#ggsave("approximation/sim_large7_biaz_K1_k0_5.png", width = 5, height = 5)
#ggsave("approximation/sim_large7_biaz_K2_k0_5.png", width = 5, height = 5)
#ggsave("approximation/sim_large7_biaz_K5_k0_5.png", width = 5, height = 5)
#ggsave("approximation/sim_large7_biaz_K10_k0_5.png", width = 5, height = 5)
