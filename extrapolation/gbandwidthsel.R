####
##  Bandwidth selection for gaussian basis
####

## collect data with subsample statistics

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
ksub_sub <- 2500
nsub_sub <- 100
mc.reps <- 10000
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))
subsub_kref_inds <- which(kref < ksub_sub)
ssi <- subsub_kref_inds

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

repno <- 20
subfun_no_accs <- function (repno) {
  set.seed(repno)
  sigma2 <- sigma2s[repno]
  sigma2_tr <- sigma2s[repno]
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  #pmat_sub <- -pdist2(ys[1:ksub, ], mu_hats[1:ksub, ])
  rSqs <- rowSums((ys - mu_hats)^2)
  #counts <- countDistEx(mu_hats, ys, rSqs)
  #accs <- sapply(Ktarg, function(k) count_acc(counts, k))    
  #accs <- count_acc(counts, Ktarg)   
  counts_sub <- countDistEx(mu_hats[1:ksub,], ys[1:ksub,], rSqs[1:ksub])
  #accs_sub <- sapply(kref, function(k) count_acc(counts_sub, k))
  accs_sub <- count_acc(counts_sub, kref)
  accs_subsub <- matrix(NA, nsub_sub, length(ssi))
  for (ii in 1:nsub_sub) {
    subinds <- sample(ksub, ksub_sub, replace = FALSE)
    counts_subsub <- countDistEx(mu_hats[subinds,], ys[subinds,], rSqs[subinds])
    accs_subsub[ii, ] <- count_acc(counts_subsub, kref[ssi])
  }
  list(accs_sub = accs_sub, accs_subsub = accs_subsub)
}  

## load fitting data

load("approximation/sim_large5_k5_raw.RData", verbose = TRUE)
mcc <- 20
sigma2_seq <- 0.005 * 1:50
mc.reps <- nrow(accsZ)
sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))
facc <- t(sapply(1:50, function(ind) colMeans(accsZ[sigma2s==sigma2s[ind], ])))
facc_rep <- facc[match(sigma2s, sigma2_seq), ]
sqrt(colMeans((accsZ - facc_rep)^2))
(max.mu <- (qnorm(1- 1/(max(kref)^2))))

## setup constructing basis 


## actually try bandwid-sel?

t1 <- proc.time()
dats <- subfun_no_accs(20)
proc.time() - t1

matplot(kref[ssi], t(dats$accs_subsub), type = "l")

## variance of difference bandwidths on predicting ksub
