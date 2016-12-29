####
##  Methods for estimating mutual information for continuous RVs
####

library(pracma)
library(FNN)

logsumexp <- function(v) max(v) + log(sum(exp(v - max(v))))

####
##  naive kernel density estimate
####

naive_kde_mi <- function(X, Y, h, mc.mult = 10) {
  mi_kdes <- sapply(1:mc.mult, function(i) {
    sampX <- X + h * randn(nrow(X), ncol(X))
    sampY <- Y + h * randn(nrow(X), ncol(X))
    dX <- pdist2(X, sampX)
    dY <- pdist2(X, sampX)
    dXY <- dX + dY
    lseX <- apply(dX, 2, logsumexp)
    lseY <- apply(dX, 2, logsumexp)
    lseXY <- apply(dXY, 2, logsumexp)
    integrands <- log(nrow(X)) - 1/2/h^2 * (lseXY - lseX - lseY)
    integrands
  })
  mean(mi_kdes)
}

####
##  cross-validated kernel density estimate
####
cv_kde_mi <- function(X, Y, h, mc.mult = 10) {
  mi_kdes <- sapply(1:mc.mult, function(i) {
    sampX <- X + h * randn(nrow(X), ncol(X))
    sampY <- Y + h * randn(nrow(X), ncol(X))
    dX <- pdist2(X, sampX)
    dY <- pdist2(X, sampX)
    diag(dX) <- -Inf
    diag(dY) <- -Inf
    dXY <- dX + dY
    lseX <- apply(dX, 2, logsumexp)
    lseY <- apply(dX, 2, logsumexp)
    lseXY <- apply(dXY, 2, logsumexp)
    integrands <- log(nrow(X) - 1) - 1/2/h^2 * (lseXY - lseX - lseY)
    integrands
  })
  mean(mi_kdes)
}


####
##  Nearest-neighbor estimate
####

# nn_entropy <- function(X) {
#   n <- nrow(X)
#   dX <- pdist(X)
#   diag(dX) <- Inf
#   nns <- apply(dX, 1, min)
#   mean(log(n * nns)) + log(2) - 0.57721566
# }

nn_mi <- function(X, Y, ...) {
#  nn_entropy(X) + nn_entropy(Y) - nn_entropy(cbind(X, Y))
  ans <- entropy(X, ...) + entropy(Y, ...) - entropy(cbind(X, Y), ...)
  pmax(0, ans)
}

