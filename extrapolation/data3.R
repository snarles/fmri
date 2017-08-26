## RESULTS FROM SHERLOCK J = 10000 (ECOSSolver)
##             20   100   400(true)
## Logistic 0.721 0.670 0.711
## SVM      0.575 0.537 0.545
## KNN2/10  0.587 0.608 0.591
## KNN5/25  0.651 0.556 0.521
## KNN9/45  0.601 0.501 0.460
## KNN15/75 0.539 0.453 0.402
## KNN25    0.504 0.384 0.341



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

load("rakesh/converted1.rda", verbose = TRUE)
# err1$TestErr <- 1 - err1$TestErr
View(err1)
View(err_knn400)
lprobs2 <- readRDS("rakesh/converted2.rds")

####
##  Stuff for 20 -> 400
####


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


dim(pmat)
ksub <- 20

sub_us <- getYs(pmat, ncl = ksub, ny = 50)
Ys <- sub_us
acck <- mean(binmom(sub_us, ksub, ksub - 1))
(ihat <- Ihat_LI(1 - acck, ksub))

K <- 400
(ah_m <- est_moment(mle_est, K-1))
(ah_c <- cm_est_moment(cm2, K-1)) 
(ah_e <- expmix_binmom(Ys, ksub, K-1))
(ah_bm <- mean(bbinom(Ys, ksub, K-1)))
(ah_i <- 1 - lineId::piK(sqrt(2 * ihat), K))

list(ac0 = ac0, ah_m = ah_m, ah_c = ah_c, ah_bm = ah_bm, ah_i = ah_i)
