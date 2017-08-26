library(parallel)
mcc <- 100

source("approximation/gaussian_identity_finsam.R")
source("extrapolation/ku_source.R")

p <- 10
sigma2 <- 0.3
sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 2000
ksub <- 1000
mc.reps <- 100

# mus <- randn(K, p)
# ys <- mus + sqrt(sigma2) * randn(K, p)
# mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
# pmat <- -pdist2(mu_hats, ys)
# accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
# plot(accs, type = "l")

nsplines <- c(100, 200, 400)
nrows <- c(125, 250, 500, 1000)
combmat <- cbind(nsplines = rep(nsplines, each = length(nrows)), 
                 nrows = rep(nrows, length(nsplines)))

basis_vecs <- list()

kseq <- function(nr) {
  interv <- floor(ksub/nr)
  seq(interv, ksub, by = interv)
}

for (ii in 1:nrow(combmat)) {
  ns <- combmat[ii, 1]
  nr <- combmat[ii, 2]
  nm <- paste0("s",ns,"r",nr)
  knts <- seq(0, 1, length.out = ns + 2)
  knts <- rev(1 - knts[-c(1, ns + 2)]^2)  
  Xmat <- spline1_moments(knts, kseq(nr))
  Xpred <- spline1_moments(knts, K)
  basis_vecs[[nm]] <- list(Xmat = Xmat, Xpred = Xpred, kseq = kseq(nr))
}


subfun <- function (repno) {
  set.seed(repno)
  mus <- randn(K, p)
  ys <- mus + sqrt(sigma2) * randn(K, p)
  mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
  pmat <- -pdist2(mu_hats, ys)
  accs <- 1 - resample_misclassification(pmat, 1:K, 1:K)
  #all_accs[repno, ] <- accs  
  pmat_sub <- pmat[1:ksub, 1:ksub]
  accs_sub <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
  preds <- numeric(length(basis_vecs))
  for (ind in 1:length(basis_vecs)) {
    Xmat <- basis_vecs[[ind]]$Xmat
    ks <- basis_vecs[[ind]]$ks
    Xpred <- basis_vecs[[ind]]$Xpred
    bt <- nnls::nnls(Xmat, 1- accs_sub[ks])
    preds[ind] <- 1 - (Xpred %*% bt$x)
  }
  list(preds = preds, accs = accs)
}  

t1 <- proc.time()
res <- mclapply(1:mc.reps, subfun, mc.cores = mcc)
(runtime <- proc.time() - t1)

# 46s

resa <- lapply(res, `[[`, "accs")
true_accs <- colMeans(do.call(rbind, resa))
plot(true_accs, type="l", ylim = c(0,1))

(facc <- true_accs[K])
all_final_preds <- do.call(rbind, lapply(res, `[[`, "preds"))
colnames(all_final_preds) <- names(basis_vecs)
resids <- all_final_preds - facc

boxplot(resids^2)
rmses <- sqrt(colMeans(resids^2))
dispmat <- matrix(rmses, length(nrows), length(nsplines))
rownames(dispmat) <- paste0("r", nrows)
colnames(dispmat) <- paste0("s", nsplines)
dispmat
#matrix(names(rmses), length(nrows), length(nsplines))
