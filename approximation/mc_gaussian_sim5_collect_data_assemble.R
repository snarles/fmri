#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")

p <- 10
sigma2_seq <- 0.005 * 1:50
# sigma2 <- 0.07
# sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 100000 ## multiple of 1000
Ktarg <- 5000 * 1:20
ksub <- 5000 ## multiple of 250
mc.reps <- 10000
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))

# mus <- randn(K, p)
# ys <- mus + sqrt(sigma2) * randn(K, p)
# mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
# pmat <- -pdist2(mu_hats, ys)
# accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
# plot(accs, type = "l")

kseq <- function(nr, ksub) {
  interv <- floor(ksub/nr)
  seq(interv, ksub, by = interv)
}

nrow <- 500
#kz <- kseq(1000, K)
kref <- kseq(nrow, ksub)

repno <- 10
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
  #accs <- sapply(Ktarg, function(k) count_acc(counts, k))    
  accs <- count_acc(counts, Ktarg)   
  counts_sub <- countDistEx(mu_hats[1:ksub,], ys[1:ksub,], rSqs[1:ksub])
  #accs_sub <- sapply(kref, function(k) count_acc(counts_sub, k))
  accs_sub <- count_acc(counts_sub, kref)
  list(accs_sub = accs_sub, accs = accs)
}  

# res <- lapply(as.numeric(args[1]):as.numeric(args[2]), subfun)
# save(res, file = paste0(args[3], "_", args[1], "_", args[2], ".RData"))

nres <- length(res)
accs_subs <- do.call(rbind, lapply(res, `[[`, "accs_sub"))
accsZ <- do.call(rbind, lapply(res, `[[`, "accs"))
sigma2s <- sigma2s[1:nres]
save(sigma2s, accs_subs, accsZ, kref, Ktarg, file = "approximation/sim_large5_k5_raw.RData")
