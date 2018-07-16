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

fprefix <- "/data/HCP_preproc/fingerprinting/Results/DC_Matrices/"
outfile <- paste0("/data/MLcore/fmri/fingerprint/results/results_dc", i1, "_", i2, ".txt")

names <- list.files(fprefix)

get_mats <- function(ind) {
  fname <- paste0(fprefix, names[ind])
  mat <- readMat(gzfile(fname))[[1]]
  mat
}


mats1 <- get_mats(i1)
mats2 <- get_mats(i2)

fcz1 <- fc_flatten(mats1)
fcz2 <- fc_flatten(mats2)

r12 <- cor(fcz1, fcz2)

cor_ident_rate <- mean(apply(r12, 1, which.max) == 1:nrow(r12))
cor_ident_rate2 <- mean(apply(t(r12), 1, which.max) == 1:nrow(r12))

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
kl_ident_rate2 <- mean(apply(-t(klds), 1, which.max) == 1:nrow(r12))


sink(outfile)
print(cor_ident_rate)
print(kl_ident_rate)
print(cor_ident_rate2)
print(kl_ident_rate2)
sink()