####
##  Identification loss using nearest-neighbor
####

library(pracma)
library(randomForest)

fitter_ols <- function(X, Y, Xte, ...) {
  Xte %*% solve(t(X) %*% X, t(X) %*% Y)
}
fitter_rf <- function(X, Y, Xte, ...) {
  Yh <- matrix(NA, nrow(Xte), ncol(Y))
  for (i in 1:ncol(Y)) {
    y <- Y[, i]
    res <- randomForest(X, y)
    yh <- predict(res, Xte)
    Yh[, i] <- yh
  }
  Yh
}

nn_loss <- function(Yhat, Yte) {
  dm <- pdist2(Yhat, Yte)
  o <- apply(dm, 2, rank)
  sum(diag(o) != 1)/nrow(Yte)
}

id_cv_loss <- function(X, Y, k, fitter = fitter_ols, mc.reps = 20, ...) {
  losses <- numeric()
  n <- nrow(X)
  for (rep.i in 1:mc.reps) {
    inds.te <- sample(n, k)
    Xtr <- X[-inds.te, ]
    Ytr <- Y[-inds.te, ]
    Xte <- X[inds.te, ]
    Yte <- Y[inds.te, ]
    Yhat <- fitter(Xtr, Ytr, Xte, ...)
    loss <- nn_loss(Yhat, Yte)
    losses[rep.i] <- loss
  }
  mean(losses)
}

