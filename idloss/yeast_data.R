
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
View(annot)
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

