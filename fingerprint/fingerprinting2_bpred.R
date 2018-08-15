library(R.matlab)
library(pracma)
library(parallel)
mcc <- 7

fcs1 = readMat(gzfile("fingerprint/fingerprint_data/FL_Matrices/All_Sub_REST1_FL_GSR.mat.gz"))
fcs2 = readMat(gzfile("fingerprint/fingerprint_data/FL_Matrices/All_Sub_REST1_FL_GSR.mat.gz"))

dim(fcs2[[1]])

fc_flatten <- function(arr) {
  apply(arr, 1, function(a) {
    a[upper.tri(a)]
  })
}

nsubs <- nrow(fcs1[[1]])

fcz1 <- fc_flatten(fcs1[[1]])
fcz2 <- fc_flatten(fcs2[[1]])

fscores <- read.csv('fingerprint/fingerprint_data/PMAT24_A_CR_338', header = FALSE)
y <- fscores[[1]]

r12 <- cor(fcz1, fcz2)


####
## Try Kld
####


mats2 <- fcs2[[1]]
mats1 <- fcs1[[1]]


nsubs1 <- 338
nsubs2 <- 338



ijs <- apply(cbind(rep(1:nsubs1, each=nsubs2), rep(1:nsubs2, nsubs1)), 1, list)
get_kl_subroutine <- function(ij) {
  i <- ij[[1]][1]
  j <- ij[[1]][2]
  vstuff <- solve(mats1[i,,], mats2[j,,])
  vstuff2 <- solve(mats2[i,,], mats1[j,,])  
  pmin(sum(diag(vstuff)) - log(det(vstuff)), sum(diag(vstuff2)) - log(det(vstuff2)))
}

get_kl_subroutine1 <- function(ij) {
  i <- ij[[1]][1]
  j <- ij[[1]][2]
  vstuff <- solve(mats1[i,,], mats2[j,,])
  sum(diag(vstuff)) - log(det(vstuff))
}

t1 <- proc.time()
kls <- mclapply(ijs, get_kl_subroutine, mc.cores = mcc)
proc.time() - t1

klds <- t(matrix(unlist(kls), nsubs2, nsubs1))
#display_results(-klds)

saveRDS(klds, file = "fingerprint/klds_fl_rest_gsr.rds")

fl_diffmat <- abs(repmat(t(y), 338, 1) - repmat(t(t(y)), 1, 338))
klds <- klds[!is.na(y), !is.na(y)]
fl_diffmat <- fl_diffmat[!is.na(y), !is.na(y)]
cor(klds[upper.tri(klds)], fl_diffmat[upper.tri(fl_diffmat)])

y2 <- y[!is.na(y)]

kl_matches <- apply(klds, 1, function(v) order(v)[2])
plot(y2, y2[kl_matches])
cor(y2, y2[kl_matches])
