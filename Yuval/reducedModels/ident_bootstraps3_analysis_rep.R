ms <- c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250)

ihatss <- list()
nvoxs <- c(100, 200, 300, 400, 500, 600)
##nvoxs <- c(20, 40, 60, 80, 100)

##fn <- 'Yuval/ident_bootstraps3_nvox'
##fn <- 'Yuval/replicates/res_nvox'
fn <- 'Yuval/reducedModels/replicates/res_voxseed01_nvox'

for (nvox in nvoxs) {
  ff <- paste0(fn, nvox, '.rds')
  ihats <- readRDS(ff)
  boxplot(ihats)
  title(paste("nvox=", nvox))
  ihatss <- c(ihatss, list(ihats))
}

bigmat <- do.call(cbind, ihatss)

boxplot(bigmat, ylim = c(0, 9))
for (i in 0:length(nvoxs))  {
  abline(v = length(ms) * i + 0.5)
  if (i > 0) {
    text(length(ms) * i - length(ms)/2, 9, paste0("n=", nvoxs[i]))
  }
}
