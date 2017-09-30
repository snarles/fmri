library(class)
library(pracma)
library(lineId)
library(parallel)
source("extrapolation/kay_method.R")
source("par2/objective_function.R")
source("extrapolation/basis_source.R")
source("approximation/gaussian_identity_finsam2.R")

# x <- read.csv("facerec/reps.csv", header = FALSE)
# x <- as.matrix(x)
# dim(x)
# labs <- read.csv("facerec/labels.csv", header = FALSE, stringsAsFactors = FALSE)
# dim(labs)
# labs[, 2] <- sapply(labs[, 2], substring, 19)
# y <- labs[, 1]
# rownames(x) <- labs[, 2]
# save(x, y, file = "facerec/lfw.RData")
load("facerec/lfw.RData")

## get 2-repeat data
yc <- table(y)
y2s <- as.numeric(names(yc)[yc > 1])
head(y2s)
filt <- y %in% y2s
xf <- x[filt, ]; yf <- y[filt]

set.seed(0)
split_inds <- sapply(unique(yf), function(x) sample(which(yf==x), 2))
length(split_inds)
split_ind_m <- matrix(split_inds, nrow = 2)

dim(split_ind_m)

set.seed(0)
nsub <- 400
clsub <- sample(ncol(split_ind_m), nsub)
tr_inds <- split_ind_m[1, clsub]
te_inds <- split_ind_m[2, clsub]
xtr <- xf[tr_inds, ]
xte <- xf[te_inds, ]
ytr <- yf[tr_inds]
yte <- yf[te_inds]

liks <- -pdist2(xte, xtr)
accs <- 1 - resample_misclassification(liks, 1:nsub, 1:nsub)

plot(accs, type = "l", ylim = c(0, 1), xlab = "faces", ylab = "accuracy",
     main = "Subset of 400", cex.lab = 1.5)

## full acc

tr_inds <- split_ind_m[1, ]
te_inds <- split_ind_m[2, ]
xtr <- xf[tr_inds, ]
xte <- xf[te_inds, ]
ytr <- yf[tr_inds]
yte <- yf[te_inds]

liks <- -pdist2(xte, xtr)
accs_full <- 1 - resample_misclassification(liks, 1:ncol(split_ind_m), 
                                       1:ncol(split_ind_m))
## set up bases

K <- length(yte) ## 1672
Ktarg <- 1:K
ksub <- nsub
ksub2 <- ksub/2
nboot <- 20 ## number of bootstraps for CV

kseq <- function(nr, ksub) {
  interv <- floor(ksub/nr)
  seq(interv, ksub, by = interv)
}

nrow <- 100
kref <- kseq(nrow, ksub)
kde_bdwids <- list("bcv", "ucv", 0.1, 0.2, 0.3, 0.4)
(max.mu <- (qnorm(1- 1/(max(kref)^2))))
fixed.gb <- c(0.4, 0.5, 0.6, 0.7)
lsub2 <- length(kref)/2
basis_vecs <- list()

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

pmat_sub <- liks[1:ksub, 1:ksub]
rSqs <- rowSums((xtr - xte)^2)
counts <- countDistEx(xtr, xte, rSqs)
accs_f2 <- count_acc(counts, Ktarg)    
counts_sub <- countDistEx(xtr[1:ksub,], xte[1:ksub,], rSqs[1:ksub])
accs_sub <- count_acc(counts_sub, kref)

boot_accs <- matrix(NA, nboot, lsub2)
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

abs(preds[, K] - accs_f2[K])

plotcols <- c("black", "green", "blue", "red")

pdf("facerec/version2_acc_plot2.pdf")
matplot(Ktarg, t(rbind(accs_full, preds[c(5, 6, 11), ])), type = "l",
        ylim = c(0, 1), xlab = "no. faces", ylab = "accuracy",
        main = "Full set (1672)", cex.lab = 1.5, col = plotcols, lwd = 2)
abline(v = nsub, lty = 2, col = "red")
legend(500, 0.3, c("true", "kdeBCV", "kdeUCV", "regCV"), col = plotcols, 
       lwd = 2, cex = 1.5)
dev.off()
