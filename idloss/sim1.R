
library(pracma)
source("idloss/mi_source.R")
source("idloss/idLoss.R")
n <- 1000
p <- 2
sgma <- 1.5

X0 <- randn(n, p)
Y0 <- X0 + sgma * randn(n, p)
## true MI
(mi_true <- p/2 * log(1 + sgma^(-2)))

## transform both X and Y using a complicated bijection

# randomly recenters v and applies sigmoid
randwarp <- function(v, nits = 2, mult = 3) {
  for (i in 1:nits) {
    v <- (v - min(v))/(max(v) - min(v))
    v <- sigmoid(v - sample(v, 1), mult * rexp(1))
  }
  v <- (v - min(v))/(max(v) - min(v))
  v
}
xs <- seq(0, 1, 0.001)
plot(xs, randwarp(xs, 20, 1), type = "l")

bmat1 <- 0.1 * randn(p) + eye(p)
bmat2 <- 0.1 * randn(p) + eye(p)
X <- X0 %*% bmat1
Y <- Y0 %*% bmat2
X <- apply(X, 2, randwarp, nits = 1, mult = 1)
Y <- apply(Y, 2, randwarp, nits = 1, mult = 1)
pairs(cbind(X, Y))

## kde estimate
(mi_kde <- naive_kde_mi(X, Y, 2))
## cv kde estimate
(mi_cv <- cv_kde_mi(X, Y, 2))

## nn estimate
(mi_nn <- nn_mi(X, Y))

## id loss
k <- 2
(idl <- id_cv_loss(X, Y, k, mc.reps = 1000))
(mi_li <- lineId::Ihat_LI(idl, k))
(mi_np <- lineId::aba_to_mi_lower(k, 1 - idl))


