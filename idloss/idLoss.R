####
##  Identification loss using nearest-neighbor
####

library(pracma)
library(glmnet)
library(randomForest)

fitter_ols <- function(Xtr, Ytr, Xte, ...) {
  Xte %*% solve(t(Xtr) %*% Xtr, t(Xtr) %*% Ytr)
}

fitter_enet <- function(Xtr, Ytr, Xte, ...) {
  Yh <- zeros(nrow(Xte), ncol(Ytr))
  for (i in 1:ncol(Yte)) {
    res <- cv.glmnet(Xtr, Ytr[, i], ...)
    Yh[, i] <- predict(res, Xte)
  }
  Yh
}

fitter_rf <- function(Xtr, Ytr, Xte, ...) {
  Yh <- matrix(NA, nrow(Xte), ncol(Ytr))
  for (i in 1:ncol(Ytr)) {
    y <- Ytr[, i]
    res <- randomForest(Xtr, y)
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
    Xtr <- X[-inds.te, , drop = FALSE]
    Ytr <- Y[-inds.te, , drop = FALSE]
    Xte <- X[inds.te, , drop = FALSE]
    Yte <- Y[inds.te, , drop = FALSE]
    Yhat <- fitter(Xtr, Ytr, Xte, ...)
    loss <- nn_loss(Yhat, Yte)
    losses[rep.i] <- loss
  }
  mean(losses)
}

