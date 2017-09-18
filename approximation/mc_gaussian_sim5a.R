## Varying sigma
## Simplified n. rows, and added kde method

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")

mcc <- 50

p <- 10
sigma2_seq <- 0.01 * 1:20
# sigma2 <- 0.07
# sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 100000 ## multiple of 1000
Ktarg <- c(10000, 20000, 50000, 100000)
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

nsplines <- c(100, 200, 400)
nrow <- 125
#kz <- kseq(1000, K)
kref <- kseq(nrow, ksub)
kde_bdwids <- list("bcv", "ucv", 0.1, 0.2, 0.3, 0.4)

basis_vecs <- list()

for (ii in 1:length(nsplines)) {
  ns <- nsplines[ii]
  nm <- paste0("spline",ns)
  knts <- seq(0, 1, length.out = ns + 2)
  knts <- rev(1 - knts[-c(1, ns + 2)]^2)  
  Xmat <- spline1_moments(knts, kref)
  Xpred <- spline1_moments(knts, Ktarg)
  basis_vecs[[nm]] <- list(Xmat = Xmat, Xpred = Xpred)
}

(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
column_names <- c(names(basis_vecs), kde.names)

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
  accs <- sapply(Ktarg, function(k) count_acc(counts, k))    
  counts_sub <- countDistEx(mu_hats[1:ksub,], ys[1:ksub,], rSqs[1:ksub])
  accs_sub <- sapply(kref, function(k) count_acc(counts_sub, k))
  preds <- matrix(NA, length(column_names), length(Ktarg))
  rownames(preds) <- column_names
  for (ind in 1:length(column_names)) {
    if (ind <= length(basis_vecs)) {
      Xmat <- basis_vecs[[ind]]$Xmat
      Xpred <- basis_vecs[[ind]]$Xpred
      bt <- nnls::nnls(Xmat, 1- accs_sub)
      preds[ind, ] <- 1 - (Xpred %*% bt$x)
    } else {
      bw <- kde_bdwids[[ind - length(basis_vecs)]]
      preds[ind, ] <- kernel_extrap(pmat_sub, Ktarg, bw = bw)
    }
  }
  list(preds = preds, accs = accs)
}  

#tester <- subfun(10)

# mc.reps <- 600
# sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))

# t1 <- proc.time()
# res <- mclapply(1:mc.reps, subfun, mc.cores = mcc)
# (runtime2 <- proc.time() - t1)

res <- mclapply(1:50, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results1.rds")
res <- mclapply(51:100, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results2.rds")
res <- mclapply(101:150, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results3.rds")
res <- mclapply(151:200, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results4.rds")
res <- mclapply(201:250, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results5.rds")
res <- mclapply(251:300, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results6.rds")
res <- mclapply(301:350, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results7.rds")
res <- mclapply(351:400, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results8.rds")
res <- mclapply(401:450, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results9.rds")
res <- mclapply(451:500, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results10.rds")
res <- mclapply(501:550, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results11.rds")
res <- mclapply(551:600, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results12.rds")
res <- mclapply(601:650, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results13.rds")
res <- mclapply(651:700, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results14.rds")
res <- mclapply(701:750, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results15.rds")
res <- mclapply(751:800, subfun, mc.cores = mcc)
saveRDS(res, "approximation/temp_results16.rds")




res <- readRDS("approximation/temp_results1.rds")
for (i in 2:16) {
  res <- c(res, readRDS(paste0("approximation/temp_results",i,".rds")))
}
length(res)

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

ind <- 1
ind <- 2
ind <- 3
ind <- 4
for (ind in 1:length(Ktarg)) {
  all_final_preds <- do.call(rbind, lapply(res, function(v) v$preds[, ind]))
  #all_final_preds <- do.call(rbind, lapply(res, `[[`, "preds"))
  sum(all_final_preds < 0)
  all_final_preds[all_final_preds < 0] <- 0
  colnames(all_final_preds) <- column_names
  all_final_predZ[[ind]] <- all_final_preds
  resids <- all_final_preds - facc[,ind][match(sigma2s, sigma2_seq)]
  rmses <- matrix(NA, length(sigma2_seq), ncol(all_final_preds))
  for (ii in 1:length(sigma2_seq)) {
    rmses[ii, ] <- sqrt(colMeans(resids[sigma2s == sigma2_seq[ii], ]^2))
  }
  colnames(rmses) <- column_names
  rmseZ[[ind]] <- rmses
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
#ggsave("approximation/sim_large5_K100_k5_1.png", width = 6, height = 4)

save(p, K, Ktarg, ksub, sigma2_seq, true_accs, all_final_predZ, file = "approximation/sim_large5_k0.5.RData")



