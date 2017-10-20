library(class)
library(pracma)
library(lineId)
library(parallel)
source("extrapolation/kay_method.R")
source("par2/objective_function.R")
source("extrapolation/basis_source.R")
source("approximation/gaussian_identity_finsam2.R")

load("facerec/temp_setup.RData", verbose = TRUE)

## set up bases
method.inds <- match(c("kde_bcv", "kde_ucv", "r.cv.gauss"), column_names)
ksub <- 100
ksub2 <- ksub/2
lsub2 <- length(kref)/2
kseq <- function(nr, ksub) {
  interv <- floor(ksub/nr)
  seq(interv, ksub, by = interv)
}
kde_bdwids <- list("bcv", "ucv", 0.1, 0.2, 0.3, 0.4)
nrow <- 100
kref <- kseq(nrow, ksub)

basis_vecs <- lapply(basis_vecs, function(set1) {
  list(Xmat = set1$Xtarg[kref, ], Xtarg = set1$Xtarg)
})
basis_sets <- lapply(basis_sets, function(set1) {
  list(Xmat = set1$Xtarg[kref, ], Xtarg = set1$Xtarg)
})
sub_basis_sets <- lapply(basis_sets, function(set1) {
  list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(kref), , drop = FALSE])
})

repno <- 1
subfun <- function(repno) {
  set.seed(repno)
  split_inds <- sapply(unique(yf), function(x) sample(which(yf==x), 2))
  length(split_inds)
  split_ind_m <- matrix(split_inds, nrow = 2)
  clsub <- sample(ncol(split_ind_m), ksub)
  tr_inds <- split_ind_m[1, clsub]
  te_inds <- split_ind_m[2, clsub]
  xtr <- xf[tr_inds, ]
  xte <- xf[te_inds, ]
  ytr <- yf[tr_inds]
  yte <- yf[te_inds]
  
  pmat_sub <- -pdist2(xte, xtr)
  counts_sub <- countDistEx(xtr, xte, rSqs)
  accs_sub <- count_acc(counts_sub, kref)
  accs <- 1 - resample_misclassification(pmat_sub, 1:ksub, 1:ksub)
  boot_accs <- matrix(NA, nboot, lsub2)
  rSqs <- rowSums((xtr - xte)^2)
  for (ii in 1:nboot) {
    subinds <- sample(ksub, ksub2, replace = FALSE)
    counts_subsub <- countDistEx(xtr[subinds,], xte[subinds,], rSqs[subinds])
    boot_accs[ii, ] <- count_acc(counts_subsub, kref[1:lsub2])
  }
  
  ## do predictions
  
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
  list(accs_sub = accs, preds = preds)
}

library(parallel)
mc.reps <- 1000; mcc <- 20
t1 <- proc.time()
res <- mclapply(1:mc.reps, subfun, mc.cores = mcc)
(runtime <- proc.time() - t1)

plotcols <- c("black", "green", "blue", "red")

mc.ind <- 5
#pdf("facerec/version2_acc_plot2.pdf")
# matplot(Ktarg, t(rbind(accs_full, res[[mc.ind]]$preds[method.inds, ])), type = "l",
#         ylim = c(0, 1), xlab = "no. faces", ylab = "accuracy",
#         main = "Full set (1672)", cex.lab = 1.5, col = plotcols, lwd = 2)
# abline(v = ksub, lty = 2, col = "red")
# legend(500, 0.3, c("true", column_names[method.inds]), col = plotcols, 
#        lwd = 2, cex = .5)
#dev.off()

library(abind)
help(abind)
predZ <- abind(lapply(res, `[[`, "preds"), along = 3)
dim(predZ)

# nplot <- 100
# tred <- rgb(1,0,0,0.2)
# for (method.ind in method.inds) {
#   pdf(paste0("facerec/repeat_", ksub, "_", column_names[method.ind], ".pdf"))
#   plot(Ktarg, accs_full, type = "l", ylim = c(0, 1), xlab = "no. faces", ylab = "accuracy",
#        main = paste0(column_names[method.ind], " (", ksub, ")"), cex.lab = 1.5, lwd = 2)
#   matplot(Ktarg, predZ[method.ind, , 1:nplot], type = "l", col = tred, lty = 1, add= TRUE)
#   lines(Ktarg, accs_full, type = "l", col = "black", lwd = 2)
#   dev.off()
# }

resids <- predZ[,K,] - accs_full[K]
rmses <- sqrt(rowMeans(resids^2))
lala <- matrix(0, 20, 3)
for (i in 1:20)
  lala[i, ] <- sqrt(rowMeans(resids[, sample(mc.reps, mc.reps, replace = TRUE)]^2))[c(5,6,11)]
rbind(rmses[c(5,6,11)],apply(lala, 2, sd))

# K 100 
# kde_bcv     kde_ucv r.cv.gauss
# [1,] 0.0525161915 0.081718288  0.1135644
# [2,] 0.0009616969 0.001420026  0.0024238

