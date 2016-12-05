source("extrapolation/ku_source.R")

## get pmat
load("rakesh/converted1.rda", verbose = TRUE)
lprobs2 <- readRDS("rakesh/converted2.rds")

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

ac0 <- 0.601
pmat <- read.table("~/github/predict_test_error/naive/tel400_logprobs/naive_tel400_20.logprobs")


dim(pmat) # 20 1000
true_ys <- rep(1:20, each = 50)


####
##  Spline basis
####

## positive

nsplines <- 1000
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
Kmax <- 400
ks <- 2:20
MM <- spline1_moments(knts, 1:Kmax)

xmat <- MM[ks, ]
avr <- get_sub_errs(pmat, true_ys, ks)

#bt <- pinv(xmat) %*% avr
bt <- nnls::nnls(xmat, avr)$x

plot(avr); lines(xmat %*% bt)

plot(1:Kmax, MM %*% bt, type = "l", main = "predicted err")

list((MM %*% bt)[Kmax], 1-ac0)

xs <- seq(0, 1, 0.01)
plot(xs, spline1_dm(knts, xs) %*% bt, type = "l", main = "K(u)")

## nonpositive

nsplines <- 5
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
Kmax <- 400
ks <- 2:20
MM <- spline1_moments(knts, 1:Kmax)

xmat <- MM[ks, ]
avr <- get_sub_errs(pmat, true_ys, ks)

bt <- pinv(xmat) %*% avr
#plot(bt)

plot(avr); lines(xmat %*% bt)

plot(1:Kmax, MM %*% bt, type = "l", main = "predicted err")

list((MM %*% bt)[Kmax], 1-ac0)

xs <- seq(0, 1, 0.01)
plot(xs, spline1_dm(knts, xs) %*% bt, type = "l", main = "K(u)")

