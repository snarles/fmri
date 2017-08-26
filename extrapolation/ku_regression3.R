source("extrapolation/ku_source.R")

## get pmat
load("rakesh/converted1.rda", verbose = TRUE)
lprobs2 <- readRDS("rakesh/converted2.rds")

(ac0 <- 1 - err1[err1$configuration=="logistic" & err1$num_classes==400, "TestErr"])
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/logistic_992191_100.logprobs")

(ac0 <- 1 - err1[err1$configuration=="svm" & err1$num_classes==400, "TestErr"])
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/svm_057733_100.logprobs")

(ac0 <- err_knn400[err_knn400$k==40, "te"])
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/knn_s618542_k100_n10")

(ac0 <- err_knn400[err_knn400$k==100, "te"])
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/knn_s618542_k100_n25")

(ac0 <- err_knn400[err_knn400$k==180, "te"])
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/knn_s618542_k100_n45")

(ac0 <- err_knn400[err_knn400$k==300, "te"])
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/knn_s618542_k100_n75")

(ac0 <- err_knn400[err_knn400$k==500, "te"])
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/knn_s618542_k100_n125")

ac0 <- 1 - 0.0140
pmat <- read.table("~/github/predict_test_error/train_tel/sub_sub_runs/deepslim_955094_100.logprobs")

ac0 <- 0.601
pmat <- read.table("~/github/predict_test_error/naive/tel400_logprobs/naive_tel400_100.logprobs")


dim(pmat) # 100 5000
true_ys <- rep(1:100, each = 50)

ks <- 2:100
avr <- get_sub_errs(pmat, true_ys, ks)
avr
res <- lm(avr ~ log(ks))
(pred_exp <- sum(res$coefficients * c(1, log(400))))

list(avr[99], mean(1-ac0), pred_exp)

####
##  Spline basis
####

## positive

nsplines <- 1000
knts <- seq(0.9, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
Kmax <- 400
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

## polynomial

d <- 5
k <- 20
ws <- exp(0 * (1:(k-1)))
ws <- ws/sum(ws)
kmat <- get_kmat(k, d)
cvec <- get_rank_prop(pmat, true_ys)

(bt <- solve(t(kmat) %*% diag(ws) %*% kmat, t(kmat) %*% diag(ws) %*% cvec))
sum(bt)

ku <- get_vande(d = d) %*% bt


yhat <- kmat %*% bt
yhat

plot(kmat[, 2], cvec, ylim = c(0, 0.5))
lines(kmat[, 2], yhat, col = "red")
lines(seq(0, 1, 0.01), ku, col = "blue")

avrisk(k, bt)
avrisk(400, bt)
