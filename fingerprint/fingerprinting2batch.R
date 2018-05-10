#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(R.matlab)
library(pracma)
library(parallel)

i1 = as.numeric(args[1])
i2 = as.numeric(args[2])
mcc = as.numeric(args[3])

fc_flatten <- function(arr) {
  apply(arr, 1, function(a) {
    a[upper.tri(a)]
  })
}

nsubs <- 338

fprefix <- "/data/HCP_preproc/fingerprinting/Results/"
outfile <- paste0("/data/MLcore/fmri/fingerprint/results", i1, "_", i2, ".txt")

#list.files(fprefix)

remove.sub332 <- logical()
remove.sub332["All_Sub_REST1_errts.mat"] <- TRUE
remove.sub332["All_Sub_REST2_errts.mat"] <- TRUE
remove.sub332["All_Sub_EMOTION_errts.mat"] <- TRUE
remove.sub332["All_Sub_GAMBLING_errts.mat"] <- FALSE
remove.sub332["All_Sub_LANGUAGE_errts.mat"] <- TRUE
remove.sub332["All_Sub_MOTOR_errts.mat"] <- FALSE
remove.sub332["All_Sub_RELATIONAL_errts.mat"] <- TRUE
remove.sub332["All_Sub_SOCIAL_errts.mat"] <- TRUE
remove.sub332["All_Sub_WM_errts.mat"] <- TRUE

get_mats <- function(ind) {
  fname <- paste0(fprefix, names(remove.sub332)[ind])
  if (remove.sub332[ind]) {
    mat <- readMat(fname)[[1]][-332,,]
  } else {
    mat <- readMat(fname)[[1]]
  }
  mat
}


mats1 <- get_mats(i1)
mats2 <- get_mats(i2)

fcz1 <- fc_flatten(mats1)
fcz2 <- fc_flatten(mats2)

r12 <- cor(fcz1, fcz2)

cor_ident_rate <- mean(apply(r12, 1, which.max) == 1:nrow(r12))

get_kl_subroutine <- function(ij) {
  i <- ij[[1]][1]
  j <- ij[[1]][2]
  vstuff <- solve(mats1[i,,], mats2[j,,])
  vstuff2 <- solve(mats2[i,,], mats1[j,,])  
  pmin(sum(diag(vstuff)) - log(det(vstuff)), sum(diag(vstuff2)) - log(det(vstuff2)))
}

ijs <- apply(cbind(rep(1:nsubs, each=nsubs), rep(1:nsubs, nsubs)), 1, list)

kls <- mclapply(ijs, get_kl_subroutine, mc.cores = mcc)
klds <- t(matrix(unlist(kls), nsubs, nsubs))

kl_ident_rate <- mean(apply(-klds, 1, which.max) == 1:nrow(r12))

sink(outfile)
print(cor_ident_rate)
print(kl_ident_rate)
sink()