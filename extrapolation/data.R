library(lineId)
source("extrapolation/constrained_mle.R")
source("extrapolation/mle_theory.R")
source("extrapolation/moment_mle.R")

subAcc <- function(pmat, cl) {
  cl_assigned <- apply(pmat, 2, function(v) which(v==max(v))[1])
  sum(cl == cl_assigned)/length(cl)
}

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
  rankconv <- apply(pmat, 2, rank)
  Us <- list()
  for (i in 1:ncl) {
    Us[[i]] <- rankconv[i, 1:ny + (ny) * (i-1)]  
  }
  Us <- do.call(c, Us)
  Us
}


# read.table('/home/snarles/github/predict_test_error/lala.txt', header = FALSE)
pmat <- readRDS("extrapolation/cifar100.rds")
cl <- rep(1:100, each = 100)
# cl_assigned <- apply(pmat, 2, function(v) which(v==max(v))[1])
# sum(cl == cl_assigned)/length(cl)

ksub <- 30

sas <- numeric(100); sas[1] <- 1
for (i in 2:100) {
  sas[i] <- subAcc(pmat[1:i, 1:(100 * i)], cl[1:(100 * i)])
}

(ihat <- Ihat_LI(1 - mean(sas[ksub]), ksub))
sub_us <- getYs(pmat[1:ksub, 1:(100 * ksub)], ncl = ksub, ny = 100)
mle_est <-  res_mixtools(sub_us, ksub)
pseq <- seq(0.7, 1, 1/10000)
# cm <- cons_mle_est(Ys, ksub, pseq, 0.01)
Ys <- sub_us
cm2 <- momk_mle_est(Ys, ksub, pseq, lbda = 0.001, mpen = 10000)
cm2[3]
c(sum(cm2$ps^ksub * cm2$gu), mean(binmom(sub_us, ksub, ksub)))
extr <- matrix(NA, 100, 5)
extr[, 1] <- sas
extr[1:ksub, ] <- sas[1:ksub]

for (i in 2:100) {
  (fa_i <- 1 - lineId::piK(sqrt(2 * ihat), i))
  (fa_m <- est_moment(mle_est, i-1))
  (fa_c <- cm_est_moment(cm2, i-1))
  extr[i, 2:4] <- c(fa_i, fa_m, fa_c)
}

for (i in 1:ksub) {
  extr[i, 5] <- mean(binmom(sub_us, ksub, i-1))
}

colnames(extr) <- c("true", "info", "MPLE", "MCMPLE", "unbiased")

View(extr)

matplot(extr, type = "l", xlab = "no. classes", ylab = "accuracy", lty = 1:5, col = 1:5, lwd = 2)
abline(v = 30, lty = 2, col = "grey", lwd = 2)
title("Extrapolating accuracy for 100 classes from 30")
legend(70, 0.5, colnames(extr), lty = 1:5, col = 1:5, lwd = 2)
