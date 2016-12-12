####
##  Testing mutual information methods on gaussian example
####

library(pracma)
source("idloss/mi_source.R")
source("idloss/idLoss.R")
n <- 1000
p <- 5
sgma <- 1.5

X <- randn(n, p)
Y <- X + sgma * randn(n, p)
## true MI
(mi_true <- p/2 * log(1 + sgma^(-2)))

## kde estimate
(mi_kde <- naive_kde_mi(X, Y, 0.3))
## cv kde estimate
(mi_cv <- cv_kde_mi(X, Y, 0.3))

## nn estimate
(mi_nn <- nn_mi(X, Y))

## id loss
k <- 20
(idl <- id_cv_loss(X, Y, k, mc.reps = 1000))
(mi_li <- lineId::Ihat_LI(idl, k))
(mi_np <- lineId::aba_to_mi_lower(k, 1 - idl))
