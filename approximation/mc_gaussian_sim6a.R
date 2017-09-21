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
mc.reps <- 200
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
#ggsave("approximation/sim_large6_K10_k5.png", width = 6, height = 4)
#ggsave("approximation/sim_large6_K20_k5.png", width = 6, height = 4)
#ggsave("approximation/sim_large6_K50_k5.png", width = 6, height = 4)
#ggsave("approximation/sim_large6_K100_k5.png", width = 6, height = 4)

save(p, K, Ktarg, ksub, sigma2_seq, true_accs, all_final_predZ, file = "approximation/sim_large6_k0.5.RData")

sapply(rmseZ, apply, 2, max)
# [,1]       [,2]       [,3]       [,4]
# spline100 0.01462340 0.02502096 0.05646269 0.11163024
# spline200 0.01440087 0.03160492 0.11151580 0.15787608
# spline400 0.01431175 0.03493707 0.12616090 0.25905367
# kde_bcv   0.03825260 0.03021578 0.03570905 0.06508511
# kde_ucv   0.02805362 0.02026482 0.05252064 0.08684443
# kde_0.1   0.04977125 0.08826948 0.15607209 0.20890419
# kde_0.2   0.06817862 0.06596947 0.05994585 0.05292582
# kde_0.3   0.20529551 0.22370422 0.25658706 0.27679651
# kde_0.4   0.36386034 0.41588105 0.47498538 0.51286588
# par2      0.01082767 0.01309258 0.02140372 0.03177893

sapply(mse_sdZ, apply, 2, max)/sqrt(mc.reps/length(sigma2_seq)) * 1/(2 * sapply(rmseZ, apply, 2, max))
# [,1]       [,2]       [,3]       [,4]
# spline100 0.43396575 0.33079833 0.19100751 0.11223058
# spline200 0.44067169 0.26188624 0.09671095 0.07935544
# spline400 0.44341581 0.23690860 0.08548448 0.04836189
# kde_bcv   0.16589869 0.27392618 0.30201865 0.19249142
# kde_ucv   0.22621163 0.40843658 0.20534399 0.14426171
# kde_0.1   0.12750445 0.09376845 0.06910139 0.05997164
# kde_0.2   0.09307986 0.12546551 0.17990901 0.23671483
# kde_0.3   0.03091181 0.03699927 0.04203173 0.04526186
# kde_0.4   0.01744091 0.01990207 0.02270554 0.02442807
# par2      0.58609632 0.63218199 0.50387501 0.39423370