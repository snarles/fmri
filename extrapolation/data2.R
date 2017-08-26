library(lineId)
source("extrapolation/constrained_mle.R")
source("extrapolation/mle_theory.R")
source("extrapolation/moment_mle.R")
source("extrapolation/bayes_binom.R")
# source("extrapolation/psuedolikelihood.R")


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

load("rakesh/converted1.rda", verbose = TRUE)
# err1$TestErr <- 1 - err1$TestErr
View(err1)
View(err_knn400)

(ac0 <- 1 - err1[err1$configuration=="logistic" & err1$num_classes==400, "TestErr"])
pmat <- lprobs$logistic_111938.logprobs

(ac0 <- 1 - err1[err1$configuration=="flat_noise" & err1$num_classes==400, "TestErr"])
pmat <- lprobs$flat_noise_067902.logprobs

(ac0 <- 1 - err1[err1$configuration=="svm" & err1$num_classes==400, "TestErr"])
pmat <- lprobs$svm_114562.logprobs

(ac0 <- err_knn400[err_knn400$k==40, "te"])
pmat <- knnprobs$knn_415154_probs_02

(ac0 <- err_knn400[err_knn400$k==100, "te"])
pmat <- knnprobs$knn_415154_probs_05

(ac0 <- err_knn400[err_knn400$k==180, "te"])
pmat <- knnprobs$knn_415154_probs_09

(ac0 <- err_knn400[err_knn400$k==300, "te"])
pmat <- knnprobs$knn_415154_probs_15

(ac0 <- err_knn400[err_knn400$k==500, "te"])
pmat <- knnprobs$knn_415154_probs_25

ac0 <- 1 - 0.0140
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/deepslim_010855_20.logprobs")

ac0 <- 0.601
pmat <- read.table("~/github/predict_test_error/naive/tel400_logprobs/naive_tel400_20.logprobs")


dim(pmat)
ksub <- 20

sub_us <- getYs(pmat, ncl = ksub, ny = 50)

mle_est <-  res_mixtools(sub_us, ksub)
pseq <- seq(0.7, 1, 1/10000)
# cm <- cons_mle_est(Ys, ksub, pseq, 0.01)
Ys <- sub_us
cm2 <- momk_mle_est(Ys, ksub, pseq, lbda = 0.001, mpen = 10000)
# cm2[3]
# c(sum(cm2$ps^ksub * cm2$gu), mean(binmom(sub_us, ksub, ksub)))
(acck <- mean(binmom(sub_us, ksub, ksub - 1)))
(ihat <- Ihat_LI(1 - acck, ksub, upper = 20))



K <- 400
(ah_m <- est_moment(mle_est, K-1))
(ah_c <- cm_est_moment(cm2, K-1)) 
(ah_e <- expmix_binmom(Ys, ksub, K-1))
(ah_bm <- mean(bbinom(Ys, ksub, K-1)))
(ah_i <- 1 - lineId::piK(sqrt(2 * ihat), K))

list(ac0 = ac0, ah_m = ah_m, ah_c = ah_c, ah_bm = ah_bm, ah_i = ah_i)


res <- fit_pm_models(Ys, ksub, us = seq(0.999, 1, 1e-5))
sum(res$gu_mm * seq(0.999, 1, 1e-5)^(19))
plot(res$gu_mm)
sum(res$gu_mm * seq(0.999, 1, 1e-5)^(399))
