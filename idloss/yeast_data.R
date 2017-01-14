
source("idloss/idLoss.R")

###
## Reading yeast data
###

dat <- read.delim("data2/yeast_cell_cycle.txt", row.names = 1)
#View(dat)
plot(dat[, 1])
dim(dat)
table(colSums(is.na(dat)))

colnames(dat)
tseries_inds <- c(7:24, 26:49, 51:67, 69:82)
colnames(dat)[tseries_inds]


image(is.na(dat))

cdc_dat <- dat[, 26:49]
colnames(cdc_dat)

sort(sample(rownames(cdc_dat), 200))
non_na <- (rowSums(is.na(cdc_dat)) == 0)
sum(non_na)
cdc_dat <- cdc_dat[non_na, ]

summary(apply(cdc_dat, 1, var))
summary(apply(cdc_dat, 1, mean))

## columns are genes
cdc_dat <- t(as.matrix(cdc_dat))
cdc_dat_s <- scale(cdc_dat)
#View(cdc_dat_s)

dim(cdc_dat_s)

ncl <- 100
res <- kmeans(t(cdc_dat_s), ncl, iter.max = 100, nstart = 10)
table(res$cluster)
matplot(cdc_dat_s[, res$cluster == 1], type = "l")
matplot(cdc_dat_s[, res$cluster == 2], type = "l")
matplot(cdc_dat_s[, res$cluster == 3], type = "l")

## pick highest-variance member of each cluster

varis <- apply(cdc_dat, 2, var)
inds <- sapply(1:ncl, function(i) which.max(varis * (res$cluster == i)))

## reduced dataset

cdc_dat_r <- cdc_dat_s[, inds]
cor(cdc_dat_r[, 1:10], method = "spearman")

res_svd <- svd(cdc_dat_r)
plot(res_svd$d)
plot(res_svd$u[, 1])
plot(res_svd$u[, 2])


###
## Obtaining cell-cycle related genes
###

library(readxl)
annot <- read_excel("data2/CellCycle98.xls", skip = 3)
#View(annot)
class(annot)
colnames(annot)
sort(table(annot$Process), decreasing = TRUE)[1:10]

## matching the rownames to the processes

rownames(dat)[1:10]
annot$ORF[1:20]

annots <- character(nrow(dat))
m_inds <- match(rownames(dat), annot$ORF)
annots[!is.na(m_inds)] <- annot$Process[m_inds[!is.na(m_inds)]]
annots[1:100]
names(annots) <- rownames(dat)
annots[is.na(annots)] <- ""

cell_cyc <- dat[annots=="cell cycle", ]
dim(cell_cyc)
matplot(t(as.matrix(cell_cyc[, tseries_inds])), type = "l")

rowSums(is.na(cell_cyc[, tseries_inds]))
colSums(is.na(cell_cyc[, tseries_inds]))
fseries_inds <- intersect(tseries_inds, which(colSums(is.na(cell_cyc))== 0))
colnames(cell_cyc)[fseries_inds]

#image(is.na(cell_cyc[, tseries_inds]))

#annots[annots == "cell cycle"]

#matplot(apply(is.na(cell_cyc) + 0, 2, jitter, factor = 0.5), type = "l")

## cross-validation prediction


tr_inds <- sample(nrow(cell_cyc), 23)
cell_cyc_tr <- t(as.matrix(cell_cyc[tr_inds, fseries_inds]))
cell_cyc_te <- t(as.matrix(cell_cyc[-tr_inds, fseries_inds]))

#View(cell_cyc_tr)

ftr <- fitter_ols
#ftr <- fitter_enet

# k <- 3
# id_cv_loss(cell_cyc_tr, cell_cyc_te, 3, fitter_ols)
# id_cv_loss(cell_cyc_tr, cell_cyc_te[, 2, drop = FALSE], 3, fitter_ols, 
#            mc.reps = 100)


errs_cyc <- sapply(1:ncol(cell_cyc_te),
                   function(i) 
                     id_cv_loss(cell_cyc_tr, cell_cyc_te[, i, drop = FALSE], 
                                k, ftr, 
                                mc.reps = 1000))
errs_cyc


id_cv_loss(cell_cyc_tr, cell_cyc_te, k, ftr,mc.reps = 1000)

## apply to random other genes

sort(table(annots), decreasing = TRUE)[1:10]
non_cell_cyc <- dat[annots=="DNA replication", ]
non_cell_cyc <- dat[annots=="transport", ]
non_cell_cyc <- dat[annots=="cytoskeleton", ]
non_cell_cyc <- dat[annots=="chromatin structure", ]


non_cell_cyc <- dat[annots!="cell cycle", ]
dim(non_cell_cyc)
filt <- (rowSums(is.na(non_cell_cyc[, fseries_inds])) == 0)
sum(filt)
ncc <- t(as.matrix(non_cell_cyc[filt, fseries_inds]))

id_cv_loss(cell_cyc_tr, ncc, k, ftr,mc.reps = 1000)


(inds_samp <- sample(ncol(ncc), 5))
errs_noncyc <- sapply(inds_samp,
                   function(i) 
                     id_cv_loss(cell_cyc_tr, ncc[, i, drop = FALSE], 
                                k, ftr, 
                                mc.reps = 100))
errs_noncyc

####
##  do the full comparisons
####

library(lineId)
cats <- c("cell cycle", "DNA replication", "transport", "cytoskeleton", "chromatin structure")
dsets <- list()
for (cat in cats) {
  non_cell_cyc <- dat[annots==cat, ]
  filt <- (rowSums(is.na(non_cell_cyc[, fseries_inds])) == 0)
  sum(filt)
  ncc <- t(as.matrix(non_cell_cyc[filt, fseries_inds]))
  dsets[[cat]] <- ncc
}

# for (ii in 1:2) {
#   non_cell_cyc <- dat[sample(nrow(dat), 20), ]
#   filt <- (rowSums(is.na(non_cell_cyc[, fseries_inds])) == 0)
#   sum(filt)
#   ncc <- t(as.matrix(non_cell_cyc[filt, fseries_inds]))
#   dsets[[paste0("noise", ii)]] <- ncc
# }



ccs <- matrix(0, length(dsets), length(dsets))
for (i in 1:(length(dsets) - 1)) {
  for (j in (i+1):length(dsets)) {
    ccs[i, j] <- cancor(dsets[[i]], dsets[[j]])$cor[1]
  }
}
ccs
# [,1] [,2] [,3]      [,4]      [,5]
# [1,]    0    1    1 1.0000000 1.0000000
# [2,]    0    0    1 0.9988279 0.9994315
# [3,]    0    0    0 0.9911839 0.9897579
# [4,]    0    0    0 0.0000000 0.9891288
# [5,]    0    0    0 0.0000000 0.0000000

library(PMA)

res <- CCA(dsets[[1]], dsets[[2]])
names(res)
res$cors

help(CCA)

pccs <- matrix(0, length(dsets), length(dsets))
pccs2 <- matrix(0, length(dsets), length(dsets))
for (i in 1:(length(dsets) - 1)) {
  for (j in (i+1):length(dsets)) {
    #pccs[i, j] <- CCA(dsets[[i]], dsets[[j]])$cors
    res_cca <- CCA.permute(dsets[[i]], dsets[[j]])
    ccor <- res_cca$cors[res_cca$penaltyxs == res_cca$bestpenaltyx]
    pccs[i, j] <- ccor
    res_cca <- CCA.permute(exp(3 * dsets[[i]]), exp(3 * dsets[[j]]))
    ccor <- res_cca$cors[res_cca$penaltyxs == res_cca$bestpenaltyx]
    pccs2[i, j] <- ccor
  }
}
pccs
pccs2

# > pccs
# [,1]      [,2]      [,3]      [,4]      [,5]
# [1,]    0 0.9570214 0.8697984 0.9223172 0.9358656
# [2,]    0 0.0000000 0.8338793 0.8856461 0.9500084
# [3,]    0 0.0000000 0.0000000 0.8376016 0.7810961
# [4,]    0 0.0000000 0.0000000 0.0000000 0.9014151
# [5,]    0 0.0000000 0.0000000 0.0000000 0.0000000
# > pccs2
# [,1]      [,2]      [,3]      [,4]      [,5]
# [1,]    0 0.9919508 0.9241681 0.8840183 0.9657301
# [2,]    0 0.0000000 0.8224505 0.6722902 0.9879212
# [3,]    0 0.0000000 0.0000000 0.5269828 0.6627958
# [4,]    0 0.0000000 0.0000000 0.0000000 0.8201928
# [5,]    0 0.0000000 0.0000000 0.0000000 0.0000000

k <- 5
ftr <- fitter_rf
infos <- matrix(0, length(dsets), length(dsets))
infos2 <- matrix(0, length(dsets), length(dsets))
mcr <- 10
for (i in 1:(length(dsets) - 1)) {
  for (j in (i+1):length(dsets)) {
    ep <- id_cv_loss(dsets[[i]], dsets[[j]], k, ftr,mc.reps = mcr)
    infos[i, j] <- lineId::aba_to_mi_lower(k, 1-ep)
    ep <- id_cv_loss(dsets[[j]], dsets[[i]], k, ftr,mc.reps = mcr)
    infos[j, i] <- lineId::aba_to_mi_lower(k, 1-ep)
    ep <- id_cv_loss(exp(3 * dsets[[i]]), exp(3 * dsets[[j]]), k, ftr,mc.reps = mcr)
    infos2[i, j] <- lineId::aba_to_mi_lower(k, 1-ep)
    ep <- id_cv_loss(exp(3 * dsets[[j]]), exp(3 * dsets[[i]]), k, ftr,mc.reps = mcr)
    infos2[j, i] <- lineId::aba_to_mi_lower(k, 1-ep)
  }
}
#infos
s_cors <- sqrt(1 - exp(-2 * infos))
best_cors <- pmax(s_cors, t(s_cors))
best_cors
s_cors2 <- sqrt(1 - exp(-2 * infos2))
best_cors2 <- pmax(s_cors2, t(s_cors2))
best_cors2

# [,1]      [,2]      [,3]      [,4]      [,5]
# [1,] 0.0000000 0.9208131 0.7884693 0.9734130 0.8068120
# [2,] 0.9208131 0.0000000 0.8525753 0.9004876 0.9213302
# [3,] 0.7884693 0.8525753 0.0000000 0.7305535 0.7105705
# [4,] 0.9734130 0.9004876 0.7305535 0.0000000 0.9344747
# [5,] 0.8068120 0.9213302 0.7105705 0.9344747 0.0000000



####
##  Get null distributions for CCA, pCCA, info
####

nullsets <- list()
for (ii in 1:100) {
  non_cell_cyc <- dat[sample(nrow(dat), 30), ]
  filt <- (rowSums(is.na(non_cell_cyc[, fseries_inds])) == 0)
  sum(filt)
  ncc <- t(as.matrix(non_cell_cyc[filt, fseries_inds]))
  nullsets[[paste0("noise", ii)]] <- ncc
}

null_ccs <- numeric()
null_pccs <- numeric()
null_infos <- numeric()

niters <- 100
k <- 3
ftr <- fitter_ols
for (jj in 1:niters) {
  d1 <- nullsets[[sample(100, 1)]]
  d2 <- nullsets[[sample(100, 1)]]
  null_ccs[jj] <- cancor(d1, d2)$cor[1]
  res_cca <- CCA.permute(dsets[[i]], dsets[[j]])
  ccor <- res_cca$cors[res_cca$penaltyxs == res_cca$bestpenaltyx]
  null_pccs[jj] <- ccor
  ep <- id_cv_loss(d1, d2, k, ftr,mc.reps = 1000)
  i1 <- lineId::aba_to_mi_lower(k, 1-ep)
  ep <- id_cv_loss(d2, d1, k, ftr,mc.reps = 1000)
  i2 <- lineId::aba_to_mi_lower(k, 1-ep)
  i3 <- pmax(i1, i2)
  ss <- sqrt(1 - exp(-2 * i3))
  null_infos[jj] <- ss
}

sapply(ccs[upper.tri(ccs)], function(v) mean(v < null_ccs))
sapply(pccs[upper.tri(ccs)], function(v) mean(v < null_pccs))
sapply(s_cors[upper.tri(ccs)], function(v) mean(v < null_infos))
