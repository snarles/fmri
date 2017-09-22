## Varying sigma
## Simplified n. rows, and added kde method

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/ku_source.R")
source("extrapolation/kay_method.R")
source("par2/objective_function.R")


mcc <- 40

p <- 10
sigma2_seq <- 0.02 * 1:20
# sigma2 <- 0.07
# sigma2_tr <- sigma2 ## equivalent to 1 nn
K <- 10000 ## multiple of 1000
Ktarg <- c(1000, 2000, 5000, 10000)
ksub <- 500 ## multiple of 250
mc.reps <- 4000
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
# spline100 0.03641128 0.07412600 0.21687262 0.27954394
# spline200 0.03647592 0.07586721 0.23229983 0.30405971
# spline400 0.03649262 0.07634008 0.23593262 0.31057313
# kde_bcv   0.09072399 0.08528164 0.08082117 0.07580294
# kde_ucv   0.06661882 0.05806094 0.04965848 0.04393042
# kde_0.1   0.07940925 0.12693789 0.20611621 0.26423671
# kde_0.2   0.02892534 0.04950055 0.09208124 0.12876894
# kde_0.3   0.09766652 0.10016040 0.10506089 0.10759312
# kde_0.4   0.20350515 0.22774690 0.26333504 0.28844270
# par2      0.03074391 0.03789307 0.04962707 0.09360574

sapply(mse_sdZ, apply, 2, max)/sqrt(mc.reps/length(sigma2_seq)) * 1/(2 * sapply(rmseZ, apply, 2, max))
# [,1]        [,2]        [,3]        [,4]
# spline100 0.011913266 0.007253754 0.006765055 0.008712828
# spline200 0.011892156 0.007087274 0.006315783 0.008010329
# spline400 0.011886714 0.007043374 0.006218535 0.007842334
# kde_bcv   0.004781286 0.006304894 0.018153105 0.032130923
# kde_ucv   0.006511333 0.009260816 0.029544911 0.055442633
# kde_0.1   0.005462554 0.004235865 0.007118097 0.009217562
# kde_0.2   0.014996446 0.010862338 0.015933269 0.018914642
# kde_0.3   0.004441412 0.005368306 0.013964809 0.022637305
# kde_0.4   0.002131530 0.002360918 0.005571440 0.008444028
# par2      0.014109374 0.014189711 0.029563610 0.026019967