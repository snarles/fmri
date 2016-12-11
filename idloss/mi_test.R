####
##  Testing mutual information methods on gaussian example
####

library(pracma)
source("idloss/mi_source.R")
n <- 200
p <- 5
sgma <- 0.5

X <- randn(n, p)
Y <- X + sgma * randn(n, p)
## true MI
(mi_true <- p/2 * log(1 + sgma^(-2)))

## kde estimate
(mi_kde <- naive_kde_mi(X, Y, 0.3))
## cv kde estimate
(mi_cv <- cv_kde_mi(X, Y, 0.3))
