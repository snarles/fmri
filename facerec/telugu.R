library(class)
library(pracma)
library(lineId)
library(parallel)
source("extrapolation/kay_method.R")
source("par2/objective_function.R")
source("extrapolation/basis_source2.R")
source("approximation/gaussian_identity_finsam2.R")

load("rakesh/converted3.rda", verbose = TRUE)
set.seed(0)

accZ20 <- lapply(lp20, function(v) 1- resample_misclassification(v, rep(1:20, each = 50), 1:20))
accZ100 <- lapply(lp100, function(v) 1- resample_misclassification(v, rep(1:100, each = 50), 1:100))

matplot(1:20, do.call(cbind, accZ20), type = "l", ylim = c(0,1))
matplot(1:100, do.call(cbind, accZ100), type = "l", ylim = c(0,1))

## predict 100 -> 400

K <- 400
ksub <- 100
ksub2 <- ksub/2
nboot <- 20 ## number of bootstraps for CV
i_chosen <- rep(1:100, each = 50)

kde_bdwids <- list("bcv", "ucv")
(max.mu <- (qnorm(1- 1/(ksub^2))))
lsub2 <- ksub2

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_basis_mat(max.mu, bd, 1:ksub, K))


(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(kde.names, "r.cv.gauss"))

subfun <- function(v) {
  res <- c(kernel_extrap2(v, i_chosen, K, bw = "ucv"),
    kernel_extrap2(v, i_chosen, K, bw = "bcv"),
    cv_reg(v, ksub2, nboot, i_chosen, basis_sets))
  names(res) <- column_names
  res
}

(preds <- lapply(lp100, subfun))
(resids <- lapply(methods, function(m) preds[[m]] - accs400[m]))
do.call(rbind, resids)
#           kde_bcv     kde_ucv   r.cv.gauss
# [1,] -0.094956221 -0.12352541 -0.002304109
# [2,] -0.001851355 -0.03314655  0.010733172
# [3,] -0.108299464 -0.19241302  0.051671117

(preds100 <- cbind(do.call(rbind, preds), accs400))
#            kde_bcv   kde_ucv r.cv.gauss accs400
# deepslim 0.8910438 0.8624746  0.9836959  0.9860
# logistic 0.7088486 0.6775535  0.7214332  0.7107
# svm      0.4369005 0.3527870  0.5968711  0.5452

## predict 20 -> 400

K <- 400
ksub <- 20
ksub2 <- ksub/2
nboot <- 20 ## number of bootstraps for CV
i_chosen <- rep(1:20, each = 50)

kde_bdwids <- list("bcv", "ucv")
(max.mu <- (qnorm(1- 1/(ksub^2))))
lsub2 <- ksub2

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_basis_mat(max.mu, bd, 1:ksub, K))


(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(kde.names, "r.cv.gauss"))

subfun <- function(v) {
  res <- c(kernel_extrap2(v, i_chosen, K, bw = "ucv"),
           kernel_extrap2(v, i_chosen, K, bw = "bcv"),
           cv_reg(v, ksub2, nboot, i_chosen, basis_sets))
  names(res) <- column_names
  res
}

(preds <- lapply(lp20, subfun))
(resids <- lapply(methods, function(m) preds[[m]] - accs400[m]))
do.call(rbind, resids)
#          kde_bcv      kde_ucv  r.cv.gauss
# [1,] -0.49572055 -0.599681566 -0.02457179
# [2,]  0.03597772 -0.009158706  0.17165853
# [3,] -0.02888129 -0.138234058  0.12726166

(preds20 <- cbind(do.call(rbind, preds), accs400))
#            kde_bcv   kde_ucv r.cv.gauss accs400
# deepslim 0.4902794 0.3863184  0.9614282  0.9860
# logistic 0.7466777 0.7015413  0.8823585  0.7107
# svm      0.5163187 0.4069659  0.6724617  0.5452

## predict 20 -> 100

K <- 100
ksub <- 20
ksub2 <- ksub/2
nboot <- 20 ## number of bootstraps for CV
i_chosen <- rep(1:20, each = 50)

kde_bdwids <- list("bcv", "ucv")
(max.mu <- (qnorm(1- 1/(ksub^2))))
lsub2 <- ksub2

(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_basis_mat(max.mu, bd, 1:ksub, K))


(kde.names <- sapply(kde_bdwids, function(v) paste0("kde_", v)))
(column_names <- c(kde.names, "r.cv.gauss"))

subfun <- function(v) {
  res <- c(kernel_extrap2(v, i_chosen, K, bw = "ucv"),
           kernel_extrap2(v, i_chosen, K, bw = "bcv"),
           cv_reg(v, ksub2, nboot, i_chosen, basis_sets))
  names(res) <- column_names
  res
}

(preds <- lapply(lp20, subfun))
(resids <- lapply(methods, function(m) preds[[m]] - accs100[m]))
do.call(rbind, resids)
#           kde_bcv     kde_ucv    r.cv.gauss
# [1,] -0.276940767 -0.34009654 -0.0002775256
# [2,] -0.007515977 -0.03287432  0.0489710983
# [3,] -0.103793049 -0.18108809  0.0610080879

(preds20_100 <- cbind(do.call(rbind, preds), accs100))
#            kde_bcv   kde_ucv r.cv.gauss accs100
# deepslim 0.7138592 0.6507035  0.9905225  0.9908
# logistic 0.8414840 0.8161257  0.8979711  0.8490
# svm      0.6544070 0.5771119  0.8192081  0.7582