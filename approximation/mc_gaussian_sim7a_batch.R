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
    geom_line() + coord_cartesian(xlim = c(0, 1)) + 
    ggtitle(paste0("Predicting K=", Ktarg[ind], " from k=", ksub))
  #ggsave("approximation/sim_large7_K10_k5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K20_k5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K50_k5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K100_k5.png", width = 6, height = 4)
  
  save(p, K, Ktarg, ksub, sigma2_seq, true_accs, all_final_predZ, file = "approximation/sim_large7_k5.RData")
  
  sapply(rmseZ, apply, 2, max)
  # [,1]       [,2]       [,3]       [,4]
  # spline100 0.011056534 0.02236130 0.04757882 0.10059901
  # spline200 0.011340475 0.02988464 0.09843228 0.16717138
  # spline400 0.011394711 0.03206902 0.12611819 0.23045890
  # kde_bcv   0.038068433 0.02861389 0.03551000 0.06466055
  # kde_ucv   0.028394427 0.01919864 0.05293413 0.08618046
  # kde_0.1   0.048618228 0.08897262 0.15416965 0.20777861
  # kde_0.2   0.068476801 0.06424987 0.05824453 0.05356398
  # kde_0.3   0.206113643 0.22683697 0.25426093 0.27443654
  # kde_0.4   0.369232181 0.41334571 0.47297803 0.51385121
  # par2      0.009199407 0.01175048 0.01617782 0.02122277
  
  sapply(mse_sdZ, apply, 2, max)/sqrt(mc.reps/length(sigma2_seq)) * 1/(2 * sapply(rmseZ, apply, 2, max))
  # [,1]        [,2]        [,3]        [,4]
  # spline100 0.083309413 0.051471601 0.031613357 0.017604320
  # spline200 0.081223524 0.038513819 0.015280824 0.010593782
  # spline400 0.080836923 0.035890458 0.011926324 0.007684568
  # kde_bcv   0.024196251 0.040224231 0.042357825 0.027388833
  # kde_ucv   0.032439935 0.059950680 0.028415057 0.020549636
  # kde_0.1   0.018945844 0.012936247 0.009756307 0.008523385
  # kde_0.2   0.013451466 0.017913993 0.025824337 0.033062836
  # kde_0.3   0.004468959 0.005074004 0.005915680 0.006453139
  # kde_0.4   0.002494672 0.002784526 0.003180119 0.003446478
  # par2      0.100127472 0.097951050 0.092974588 0.083447011
}