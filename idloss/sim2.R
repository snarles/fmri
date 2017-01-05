## linear simulation


library(pracma)
source("idloss/mi_source.R")
source("idloss/idLoss.R")
n <- 1000
p <- 2
sgma <- 0.1

d_ex <- 0

set.seed(0)

X0 <- randn(n, p)
Y0 <- X0 + sgma * randn(n, p)
## true MI
(mi_true <- p/2 * log(1 + sgma^(-2)))

bmat1 <- 0.1 * randn(p) + eye(p)
bmat2 <- 0.1 * randn(p) + eye(p)
X <- X0 %*% bmat1
Y <- Y0 %*% bmat2
colnames(X) <- paste0("X", 1:p)
colnames(Y) <- paste0("Y", 1:p)
pairs(cbind(X, Y))

if (d_ex > 0) {
  X <- cbind(X, randn(n, d_ex))
  #Y <- cbind(Y, randn(n, d_ex))
}
colnames(X) <- paste0("X", 1:ncol(X))



## nn estimate
(mi_nn <- nn_mi(X, Y)[1])

## id loss using OLS
k <- 10
(idl <- id_cv_loss(X, Y, k, mc.reps = 1000))
(mi_li <- lineId::Ihat_LI(idl, k))
(mi_np_ols_10 <- lineId::aba_to_mi_lower(k, 1 - idl))

k <- 20
(idl <- id_cv_loss(X, Y, k, mc.reps = 1000))
(mi_li <- lineId::Ihat_LI(idl, k))
(mi_np_ols_20 <- lineId::aba_to_mi_lower(k, 1 - idl))

## id loss using enet
k <- 10
(idl <- id_cv_loss(X, Y, k, mc.reps = 1000, fitter = fitter_enet, alpha = 0.5))
(mi_li <- lineId::Ihat_LI(idl, k))
(mi_np_enet_10 <- lineId::aba_to_mi_lower(k, 1 - idl))

k <- 20
(idl <- id_cv_loss(X, Y, k, mc.reps = 1000, fitter = fitter_enet, alpha = 0.5))
(mi_li <- lineId::Ihat_LI(idl, k))
(mi_np_enet_20 <- lineId::aba_to_mi_lower(k, 1 - idl))



## id loss using RF
k <- 10
(idl <- id_cv_loss(X, Y, k, mc.reps = 50, fitter = fitter_rf))
(mi_li <- lineId::Ihat_LI(idl, k))
(mi_np_rf_10 <- lineId::aba_to_mi_lower(k, 1 - idl))

k <- 20
(idl <- id_cv_loss(X, Y, k, mc.reps = 50, fitter = fitter_rf))
(mi_np_rf_20 <- lineId::aba_to_mi_lower(k, 1 - idl))

list(mi_true = mi_true, mi_nn = mi_nn, 
     mi_np_ols_10 = mi_np_ols_10, mi_np_ols_20 = mi_np_ols_20,
     mi_np_enet_10 = mi_np_enet_10, mi_np_enet_20 = mi_np_enet_20,
     mi_np_rf_10 = mi_np_rf_10, mi_np_rf_20 = mi_np_rf_20)
