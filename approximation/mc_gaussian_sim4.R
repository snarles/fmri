## Varying sigma

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")

mcc <- 200

p <- 10
sigma2_seq <- 0.02 * 1:20
# sigma2 <- 0.07
# sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 10000 ## multiple of 1000
ksub <- 5000 ## multiple of 250
mc.reps <- 800
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

nsplines <- c(400, 600, 800)
nrows <- c(125, 250)
kz <- kseq(1000, K)
kref <- kseq(250, ksub)

combmat <- cbind(nsplines = rep(nsplines, each = length(nrows)), 
                 nrows = rep(nrows, length(nsplines)))

basis_vecs <- list()


for (ii in 1:nrow(combmat)) {
  ns <- combmat[ii, 1]
  nr <- combmat[ii, 2]
  nm <- paste0("s",ns,"r",nr)
  knts <- seq(0, 1, length.out = ns + 2)
  knts <- rev(1 - knts[-c(1, ns + 2)]^2)  
  ks <- kseq(nr, ksub)
  Xmat <- spline1_moments(knts, ks)
  Xpred <- spline1_moments(knts, K)
  basis_vecs[[nm]] <- list(Xmat = Xmat, Xpred = Xpred, kseq = ks)
}

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
  accs <- sapply(kz, function(k) count_acc(counts, k))    
  counts_sub <- countDistEx(mu_hats[1:ksub,], ys[1:ksub,], rSqs[1:ksub])
  accs_sub <- sapply(kref, function(k) count_acc(counts_sub, k))
  preds <- numeric(length(basis_vecs) + 1)
  for (ind in 1:length(basis_vecs)) {
    Xmat <- basis_vecs[[ind]]$Xmat
    ks <- basis_vecs[[ind]]$ks
    Xpred <- basis_vecs[[ind]]$Xpred
    bt <- nnls::nnls(Xmat, 1- accs_sub[match(ks, kref)])
    preds[ind] <- 1 - (Xpred %*% bt$x)
  }
  ## kay et al method
  preds[length(preds)] <- kernel_extrap(pmat_sub, K)
  list(preds = preds, accs = accs)
}  

t1 <- proc.time()
res <- mclapply(1:mc.reps, subfun, mc.cores = mcc)
(runtime2 <- proc.time() - t1)

resa <- do.call(rbind, lapply(res, `[[`, "accs"))
true_accs <- matrix(NA, length(sigma2_seq), ncol(resa))
for (ii in 1:length(sigma2_seq)) {
  true_accs[ii, ] <- colMeans(resa[sigma2s == sigma2_seq[ii], ])
}
matplot(kz, t(true_accs), type = "l", ylim = c(0, 1))


(facc <- true_accs[, match(K, kz)])
all_final_preds <- do.call(rbind, lapply(res, `[[`, "preds"))
sum(all_final_preds < 0)
all_final_preds[all_final_preds < 0] <- 0
colnames(all_final_preds) <- c(names(basis_vecs), "kde")
resids <- all_final_preds - facc[match(sigma2s, sigma2_seq)]

rmses <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
for (ii in 1:length(sigma2_seq)) {
  rmses[ii, ] <- sqrt(colMeans(resids[sigma2s == sigma2_seq[ii], ]^2))
}
colnames(rmses) <- c(names(basis_vecs), "kde")


resids2 <- (all_final_preds > .5) - (facc[match(sigma2s, sigma2_seq)] > .5)

acc_class <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
for (ii in 1:length(sigma2_seq)) {
  acc_class[ii, ] <- colMeans(abs(resids[sigma2s == sigma2_seq[ii], ]))
}


# temp <- data.frame(sigma2_seq, rmses)
# temp2 <- melt(data = temp, id.vars = "sigma2_seq")
# colnames(temp2)[3] <- "rmse"
# library(ggplot2)
# ggplot(data = temp2, aes(x = sigma2_seq, y = rmse, colour = variable)) +
#   geom_line() + ggtitle("Predicting K=100,000 from k=25,000")
# ggsave("approximation/simulation_large_01.png", width = 6, height = 4)


temp <- data.frame(true_acc = true_accs[, ncol(true_accs)], rmses)
library(reshape2)
temp2 <- melt(data = temp, id.vars = "true_acc")
colnames(temp2)[3] <- "rmse"

library(ggplot2)
ggplot(data = temp2, aes(x = true_acc, y = rmse, colour = variable)) +
  geom_line() + coord_cartesian(xlim = c(0, 1)) + 
  ggtitle("Predicting K=10,000 from k=5,000")
ggsave("approximation/sim_large4_K10_k5_1.png", width = 6, height = 4)

save(p, K, ksub, sigma2_seq, true_accs, all_final_preds, file = "approximation/sim_large4_K10_k5.RData")



temp <- data.frame(true_acc = true_accs[, ncol(true_accs)], rmses[, 1:6])
library(reshape2)
temp2 <- melt(data = temp, id.vars = "true_acc")
colnames(temp2)[3] <- "rmse"

library(ggplot2)
ggplot(data = temp2, aes(x = true_acc, y = rmse, colour = variable)) +
  geom_line() + coord_cartesian(xlim = c(0, 1)) + 
  ggtitle("Predicting K=10,000 from k=5,000")
ggsave("approximation/sim_large4_K10_k5_2.png", width = 6, height = 4)

