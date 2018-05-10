library(R.matlab)
library(pracma)
library(parallel)
mcc <- 20

#fcs1 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST1_RP.mat")
#fcs2 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST2_RP.mat")

#fcs1 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST1_errts.mat")
#fcs2 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST2_errts.mat")

fcs1 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST1_errts.mat")
fcs2 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST1_dropout_errts.mat")

#subject 332 for motor

dim(fcs2[[1]])

#fcs1 = readMat("~/Desktop/Results/All_Sub_REST1_TP.mat")
#fcs2 = readMat("~/Desktop/Results/All_Sub_REST2_TP.mat")

#fcs1 = readMat("~/Desktop/Results/All_Sub_Rest1.mat")
#fcs2 = readMat("~/Desktop/Results/All_Sub_Rest2.mat")

#fc = fcs1[[1]][1,,]
#hist(fc[upper.tri(fc)])
#upper.tri(matrix(1:4, 2, 2))

fc_flatten <- function(arr) {
  apply(arr, 1, function(a) {
    a[upper.tri(a)]
  })
}

nsubs <- nrow(fcs1[[1]])
## shuffle subjects
shufind <- sample(nsubs, nsubs, replace = FALSE)

fcz1 <- fc_flatten(fcs1[[1]])
fcz2 <- fc_flatten(fcs2[[1]])

#r12 <- atanh(cor(fcz1, fcz2))

r12 <- cor(fcz1, fcz2)


####
## Try Kld
####


#mats1 <- fcs1[[1]]
mats2 <- fcs2[[1]]
mats1 <- fcs1[[1]][-332,,]


#nsubs1 <- dim(mats1)[1]
#nsubs1
#nsubs2 <- dim(mats2)[1]
#nsubs2

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

image(klds)
image(klds[200:300,200:300])
image(klds[300:338,300:338])
#image(apply(klds, 1, rank))
