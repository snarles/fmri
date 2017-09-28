#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

## Varying sigma
## Simplified n. rows, and added kde method and par2

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")
source("par2/objective_function.R")
source("extrapolation/basis_source.R")

p <- 10
sigma2_seq <- 0.005 * 1:50
K <- 100000 ## multiple of 1000
Ktarg <- c(10000, 20000, 50000, 100000)
ksub <- 5000 ## multiple of 250
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
  basis_vecs[[nm]] <- get_basis_mat(max.mu, kernel_sigma = gb)
}

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_basis_mat(max.mu, bd))

sub_basis_sets <- lapply(basis_sets, function(set1) {
  list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(kref), , drop = FALSE])
})

(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(names(basis_vecs), kde.names, "r.cv.gauss", "par2"))

repno <- 25
subfun <- function (repno) {
  set.seed(repno)
  sigma2 <- sigma2s[repno]
  sigma2_tr <- sigma2s[repno]
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  pmat_sub <- -pdist2(ys[1:ksub, ], mu_hats[1:ksub, ])
  rSqs <- rowSums((ys - mu_hats)^2)
  counts <- countDistEx(mu_hats, ys, rSqs)
  accs <- count_acc(counts, Ktarg)    
  counts_sub <- countDistEx(mu_hats[1:ksub,], ys[1:ksub,], rSqs[1:ksub])
  accs_sub <- count_acc(counts_sub, kref)
  boot_accs <- matrix(NA, nboot, lsub2)
  for (ii in 1:nboot) {
    subinds <- sample(ksub, ksub2, replace = FALSE)
    counts_subsub <- countDistEx(mu_hats[subinds,], ys[subinds,], rSqs[subinds])
    boot_accs[ii, ] <- count_acc(counts_subsub, kref[1:lsub2])
  }
  preds <- matrix(NA, length(column_names), length(Ktarg))
  rownames(preds) <- column_names
  for (ind in 1:length(column_names)) {
    if (ind <= length(basis_vecs)) {
      Xmat <- basis_vecs[[ind]]$Xmat
      Xpred <- basis_vecs[[ind]]$Xtarg
      bt <- nnls::nnls(Xmat, accs_sub)
      preds[ind, ] <- (Xpred %*% bt$x)
    } else if(ind <= length(basis_vecs) + length(kde_bdwids)) {
      bw <- kde_bdwids[[ind - length(basis_vecs)]]
      preds[ind, ] <- kernel_extrap(pmat_sub, Ktarg, bw = bw)
    } else if (column_names[ind] == "par2") {
      preds[ind, ] <- par2_extrapolate(kref, accs_sub, Ktarg)
    } else if (column_names[ind] == "r.cv.gauss") {
      all_sub_preds <- t(apply(boot_accs, 1, bdwid_all_preds, basis_sets = sub_basis_sets))
      cv_curve <- sqrt(colMeans((all_sub_preds - accs_sub[length(kref)])^2))
      sel_ind <- which.min(cv_curve)
      Xmat <- basis_sets[[sel_ind]]$Xmat
      Xpred <- basis_sets[[sel_ind]]$Xtarg
      bt <- nnls::nnls(Xmat, accs_sub)
      preds[ind, ] <- (Xpred %*% bt$x)
    }
  }
  list(preds = preds, accs = accs)
}  

res <- lapply(as.numeric(args[1]):as.numeric(args[2]), subfun)
save(res, file = paste0(args[3], "_", args[1], "_", args[2], ".RData"))

if (FALSE) {
  sigma2s <- sigma2s[1:length(res)]
  
  resa <- do.call(rbind, lapply(res, `[[`, "accs"))
  true_accs <- matrix(NA, length(sigma2_seq), ncol(resa))
  for (ii in 1:length(sigma2_seq)) {
    true_accs[ii, ] <- colMeans(resa[sigma2s == sigma2_seq[ii], ])
  }
  matplot(Ktarg, t(true_accs), type = "l", ylim = c(0, 1))
  
  
  (facc <- true_accs)
  
  rmseZ <- list()
  all_final_predZ <- list()
  mse_sdZ <- list()
  
  for (ind in 1:length(Ktarg)) {
    all_final_preds <- do.call(rbind, lapply(res, function(v) v$preds[, ind]))
    #all_final_preds <- do.call(rbind, lapply(res, `[[`, "preds"))
    sum(all_final_preds < 0)
    all_final_preds[all_final_preds < 0] <- 0
    colnames(all_final_preds) <- column_names
    all_final_predZ[[ind]] <- all_final_preds
    resids <- all_final_preds - facc[,ind][match(sigma2s, sigma2_seq)]
    rmses <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
    mse_sd <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
    for (ii in 1:length(sigma2_seq)) {
      rmses[ii, ] <- sqrt(colMeans(resids[sigma2s == sigma2_seq[ii], ]^2))
      mse_sd[ii, ] <- sd(resids[sigma2s == sigma2_seq[ii], ]^2)
    }
    colnames(rmses) <- column_names
    colnames(mse_sd) <- column_names
    rmseZ[[ind]] <- rmses
    mse_sdZ[[ind]] <- mse_sd
  }
  
  
  
  
  ind <- 1
  ind <- 2
  ind <- 3
  ind <- 4
  rmses <- rmseZ[[ind]]
  temp <- data.frame(true_acc = true_accs[, ind], rmses)
  library(reshape2)
  temp2 <- melt(data = temp, id.vars = "true_acc")
  colnames(temp2)[3] <- "rmse"
  
  library(ggplot2)
  ggplot(data = temp2, aes(x = true_acc, y = rmse, colour = variable)) +
    geom_line(aes(linetype=variable)) + coord_cartesian(xlim = c(0, 1)) + 
    scale_linetype_manual(values = c(1,1,1,1,1, 2,2,2,2,2,2,2, 1, 3))+
    ggtitle(paste0("Predicting K=", Ktarg[ind], " from k=", ksub))
  #ggsave("approximation/sim_large7_K10_k5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K20_k5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K50_k5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K100_k5.png", width = 6, height = 4)
  
  save(p, K, Ktarg, ksub, sigma2_seq, true_accs, rmseZ, file = "approximation/sim_large7_k5.RData")
  
  sapply(rmseZ, apply, 2, max)
#   [,1]       [,2]       [,3]       [,4]
#   r.gauss0.4  0.010574810 0.01630578 0.03191626 0.04708675
#   r.gauss0.45 0.010466451 0.01473328 0.02763747 0.04299542
#   r.gauss0.5  0.010346165 0.01446117 0.02577095 0.03689831
#   r.gauss0.55 0.010360535 0.01440899 0.02274998 0.03038103
#   r.gauss0.6  0.009373865 0.01526569 0.02451539 0.03787950
#   kde_bcv     0.038270049 0.02841774 0.03576282 0.06475881
#   kde_ucv     0.028146638 0.01909893 0.05272348 0.08606347
#   kde_0.1     0.049663395 0.08827276 0.15375912 0.20738325
#   kde_0.15    0.018357575 0.03684732 0.07908119 0.11544856
#   kde_0.2     0.068497574 0.06474735 0.05779945 0.05359504
#   kde_0.25    0.133080165 0.14040762 0.14958683 0.15447090
#   kde_0.3     0.206398817 0.22680190 0.25403324 0.27408834
#   r.cv.gauss  0.009930697 0.01529758 0.03246833 0.05474862
#   par2        0.009238257 0.01134106 0.01670042 0.02309860
  
  sapply(mse_sdZ, apply, 2, max)/sqrt(mc.reps/length(sigma2_seq)) * 1/(2 * sapply(rmseZ, apply, 2, max))
# [,1]        [,2]         [,3]         [,4]
# r.gauss0.4  0.0120795284 0.009377828 0.0059458570 0.0046402538
# r.gauss0.45 0.0122045877 0.010378735 0.0068663868 0.0050818085
# r.gauss0.5  0.0123464801 0.010574029 0.0073637007 0.0059215302
# r.gauss0.55 0.0123293545 0.010612316 0.0083415236 0.0071918067
# r.gauss0.6  0.0136271122 0.010016766 0.0077408321 0.0057681465
# kde_bcv     0.0033378247 0.005380892 0.0053063361 0.0033739730
# kde_ucv     0.0045383294 0.008006353 0.0035993361 0.0025387598
# kde_0.1     0.0025720899 0.001732276 0.0012342001 0.0010535783
# kde_0.15    0.0069583653 0.004149902 0.0023996797 0.0018925701
# kde_0.2     0.0018648648 0.002361684 0.0032832414 0.0040767671
# kde_0.25    0.0009598629 0.001089063 0.0012686246 0.0014144702
# kde_0.3     0.0006188927 0.000674213 0.0007470264 0.0007971681
# r.cv.gauss  0.0128630166 0.009995879 0.0058447579 0.0039908677
# par2        0.0138271442 0.013483115 0.0113631618 0.0094592107
}