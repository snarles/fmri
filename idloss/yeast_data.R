
###
## Reading yeast data
###

dat <- read.delim("data2/yeast_cell_cycle.txt", row.names = 1)
#View(dat)
plot(dat[, 1])
dim(dat)
table(colSums(is.na(dat)))

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

colnames(cdc_dat)[1:10]
annot$ORF[1:20]

annots <- character(ncol(cdc_dat))
m_inds <- match(colnames(cdc_dat), annot$ORF)
annots[!is.na(m_inds)] <- annot$Process[m_inds[!is.na(m_inds)]]
annots[1:10]
names(annots) <- colnames(cdc_dat)
