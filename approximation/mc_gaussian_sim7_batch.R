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
  #ggsave("approximation/sim_large7_K1_k0.5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K2_k0.5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K5_k0.5.png", width = 6, height = 4)
  #ggsave("approximation/sim_large7_K10_k0.5.png", width = 6, height = 4)
  
  save(p, K, Ktarg, ksub, sigma2_seq, true_accs, rmseZ, file = "approximation/sim_large7_k0.5.RData")
  
  sapply(rmseZ, apply, 2, max)
#   [,1]       [,2]       [,3]       [,4]
#   r.gauss0.4  0.03103291 0.04464645 0.07998869 0.11608727
#   r.gauss0.45 0.03005329 0.04370776 0.07240579 0.09822427
#   r.gauss0.5  0.02980594 0.04005612 0.06391940 0.08963934
#   r.gauss0.55 0.02939276 0.03868845 0.05619325 0.07552113
#   r.gauss0.6  0.02865736 0.03602823 0.05101958 0.06579360
#   kde_bcv     0.09084960 0.08776910 0.07931951 0.07561510
#   kde_ucv     0.06709211 0.05909782 0.05108257 0.04544760
#   kde_0.1     0.07847689 0.12777010 0.20496379 0.26768654
#   kde_0.15    0.05381567 0.08759146 0.14817503 0.19880623
#   kde_0.2     0.02941612 0.05142659 0.09169807 0.12739544
#   kde_0.25    0.05526466 0.04882415 0.04526504 0.06416641
#   kde_0.3     0.09991109 0.10181079 0.10408667 0.10744369
#   r.cv.gauss  0.03201258 0.04398329 0.07346593 0.09768539
#   par2        0.02801916 0.03582881 0.05461481 0.08820143
  
  sapply(mse_sdZ, apply, 2, max)/sqrt(mc.reps/length(sigma2_seq)) * 1/(2 * sapply(rmseZ, apply, 2, max))
# [,1]         [,2]         [,3]         [,4]
# r.gauss0.4  0.0012481642 0.0011925503 0.0016649155 0.0019688649
# r.gauss0.45 0.0012888495 0.0012181619 0.0018392786 0.0023269215
# r.gauss0.5  0.0012995453 0.0013292135 0.0020834742 0.0025497751
# r.gauss0.55 0.0013178132 0.0013762022 0.0023699363 0.0030264397
# r.gauss0.6  0.0013516307 0.0014778169 0.0026102610 0.0034738963
# kde_bcv     0.0004263548 0.0006066273 0.0016789616 0.0030226789
# kde_ucv     0.0005773282 0.0009009322 0.0026070421 0.0050290919
# kde_0.1     0.0004935742 0.0004167104 0.0006497461 0.0008538351
# kde_0.15    0.0007197563 0.0006078576 0.0008987642 0.0011496629
# kde_0.2     0.0013167667 0.0010353230 0.0014523143 0.0017941000
# kde_0.25    0.0007008850 0.0010905082 0.0029421032 0.0035619907
# kde_0.3     0.0003876864 0.0005229616 0.0012794570 0.0021272553
# r.cv.gauss  0.0012099672 0.0012105311 0.0018127371 0.0023397579
# par2        0.0013824173 0.0014860426 0.0024384305 0.0025913429
}