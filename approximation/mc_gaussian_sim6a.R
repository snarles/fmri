## Varying sigma
## Simplified n. rows, and added kde method

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")
source("par2/objective_function.R")


mcc <- 40

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
column_names <- c(names(basis_vecs), kde.names, "par2")

repno <- 2
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
  preds <- matrix(NA, length(column_names), length(Ktarg))
  rownames(preds) <- column_names
  for (ind in 1:length(column_names)) {
    if (ind <= length(basis_vecs)) {
      Xmat <- basis_vecs[[ind]]$Xmat
      Xpred <- basis_vecs[[ind]]$Xpred
      bt <- nnls::nnls(Xmat, 1- accs_sub)
      preds[ind, ] <- 1 - (Xpred %*% bt$x)
    } else if(ind <= length(basis_vecs) + length(kde_bdwids)) {
      bw <- kde_bdwids[[ind - length(basis_vecs)]]
      preds[ind, ] <- kernel_extrap(pmat_sub, Ktarg, bw = bw)
    } else {
      preds[ind, ] <- par2_extrapolate(kref, accs_sub, Ktarg)
    }
  }
  list(preds = preds, accs = accs)
}  

#res <- lapply(1:60, subfun)

#tester <- subfun(10)

# mc.reps <- 600
# sigma2s <- rep(sigma2_seq, floor(mc.reps/length(sigma2_seq)))

t1 <- proc.time()
res <- mclapply(1:mc.reps, subfun, mc.cores = mcc)
(runtime2 <- proc.time() - t1)

sigma2s <- rep(sigma2_seq, floor(length(res)/length(sigma2_seq)))

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
#ggsave("approximation/sim_large5_K1_k0.5.png", width = 6, height = 4)
#ggsave("approximation/sim_large5_K2_k0.5.png", width = 6, height = 4)
#ggsave("approximation/sim_large5_K5_k0.5.png", width = 6, height = 4)
#ggsave("approximation/sim_large5_K10_k0.5.png", width = 6, height = 4)

save(p, K, Ktarg, ksub, sigma2_seq, true_accs, all_final_predZ, file = "approximation/sim_large6_k0.5.RData")

sapply(rmseZ, apply, 2, max)
# [,1]       [,2]       [,3]       [,4]
# spline100 0.04013241 0.08099923 0.22341812 0.28150591
# spline200 0.04021022 0.08272697 0.23597806 0.30720315
# spline400 0.04022913 0.08318922 0.23807945 0.31295652
# kde_bcv   0.09175596 0.08594437 0.07947018 0.07786111
# kde_ucv   0.06759396 0.05752887 0.05065847 0.04550220
# kde_0.1   0.07779276 0.12662877 0.20539093 0.26503722
# kde_0.2   0.02942975 0.04814230 0.09120395 0.12756712
# kde_0.3   0.09954355 0.09898549 0.10357383 0.10946501
# kde_0.4   0.20509675 0.22677891 0.26204699 0.29067162

sapply(mse_sdZ, apply, 2, max)/sqrt(mc.reps/length(sigma2_seq)) * 1/(2 * sapply(rmseZ, apply, 2, max))
# [,1]        [,2]        [,3]       [,4]
# spline100 0.016021576 0.009681700 0.009936950 0.01256630
# spline200 0.015990571 0.009479499 0.009408056 0.01151515
# spline400 0.015983056 0.009426825 0.009325016 0.01130345
# kde_bcv   0.007007550 0.009124627 0.027936199 0.04543332
# kde_ucv   0.009512455 0.013631596 0.043824753 0.07774325
# kde_0.1   0.008265351 0.006192986 0.010809118 0.01334714
# kde_0.2   0.021848115 0.016289423 0.024342091 0.02773041
# kde_0.3   0.006459328 0.007922477 0.021434901 0.03231616
# kde_0.4   0.003135030 0.003458039 0.008472125 0.01217005