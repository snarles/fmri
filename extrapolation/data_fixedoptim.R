library(lineId)
source("extrapolation/constrained_mle.R")
source("extrapolation/mle_theory.R")
source("extrapolation/moment_mle.R")
source("extrapolation/bayes_binom.R")

getUs <- function(pmat, ncl, ny) {
  rankconv <- (apply(pmat, 2, rank) - 0.5)/nrow(pmat)
  Us <- list()
  for (i in 1:ncl) {
    Us[[i]] <- rankconv[i, 1:ny + (ny) * (i-1)]  
  }
  Us <- do.call(c, Us)
  Us
}

getYs <- function(pmat, ncl, ny) {
  rankconv <- apply(pmat, 2, function(v) rank(v, ties.method = "random"))
  Us <- list()
  for (i in 1:ncl) {
    Us[[i]] <- rankconv[i, 1:ny + (ny) * (i-1)]  
  }
  Us <- do.call(c, Us)
  Us
}

pmats <- list()
pmats[[1]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/logistic_992191_100.logprobs")
pmats[[2]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/svm_057733_100.logprobs")
pmats[[3]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/knn_s618542_k100_n10")
pmats[[4]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/deepslim_955094_100.logprobs")

pmats0 <- list()
pmats0[[1]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/logistic_389795_20.logprobs")
pmats0[[2]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/svm_925163_20.logprobs")
pmats0[[3]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/knn_s618542_k020_n02")
pmats0[[4]] <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/deepslim_010855_20.logprobs")

res1 <- numeric()
res0 <- numeric()
ksub <- 20
K <- 100
for (i in 1:4) {
  pmat <- pmats[[i]]
  sub_us <- getYs(pmat, ncl = 100, ny = 50)
  res1[i] <- mean(binmom(sub_us, K, ksub - 1))
  pmat <- pmats0[[i]]
  sub_us <- getYs(pmat, ncl = 20, ny = 50)
  res0[i] <- mean(binmom(sub_us, 20, ksub - 1))
}

res1
res0

plot(res0, res1, pch = 1, xlim = c(0.88, 1), ylim = c(0.88, 1), cex = 3)
abline(0, 1, lty = 2, lwd = 2)
text(res0[1], res1[1], "L")
text(res0[2], res1[2], "S")
text(res0[3], res1[3], expression(epsilon))
text(res0[4], res1[4], "D")
