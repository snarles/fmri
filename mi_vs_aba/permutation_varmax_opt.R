

source("mi_vs_aba/permutation_varmax.R")

mono_proj <- function(v) isoreg(v)$yf

itersolve <- function(Gmat, v0 = NULL, nits = 10) {
  d <- nrow(Gmat)
  if (is.null(v0)) {
    v <- runif(d); v <- sort(v/sum(v), FALSE)
  } else {
    v <- v0
  }
  obf <- numeric()
  for (i in 1:nits) {
    v
    (obf[i] <- t(v) %*% Gmat %*% v)
    v <- Gmat %*% v
    v <- pmax(0, v)
    v <- mono_proj(v)
    if (sum(v) < 1) {
      v <- v/sum(v)
    }
    sum(v)
    while (sum(v) > 1) {
      smalv <- min(v[v > 0])
      nonz <- sum(v > 0)
      exs <- sum(v) - 1
      delt <- pmin(smalv, exs/nonz)
      v <- v - delt
      v <- pmax(0, v)
    }
  }
  list(v = v, obf = c(obf, t(v) %*% Gmat %*% v))
}

gradsolve <- function(Gmat, v0 = NULL, nits = 10, eps = 0.1) {
  d <- nrow(Gmat)
  if (is.null(v0)) {
    v <- runif(d); v <- sort(v/sum(v), FALSE)
  } else {
    v <- v0
  }
  obf <- numeric()
  for (i in 1:nits) {
    v_old <- v
    (obf[i] <- t(v) %*% Gmat %*% v)
    dv <- Gmat %*% v
    v <- v + eps * runif(1) * dv
    v <- pmax(0, v)
    v <- mono_proj(v)
    if (sum(v) < 1) {
      v <- v/sum(v)
    }
    sum(v)
    while (sum(v) > 1) {
      smalv <- min(v[v > 0])
      nonz <- sum(v > 0)
      exs <- sum(v) - 1
      delt <- pmin(smalv, exs/nonz)
      v <- v - delt
      v <- pmax(0, v)
    }
    if (t(v) %*% Gmat %*% v < obf[i]) {
      v <- v_old
    }
  }
  list(v = v, obf = c(obf, t(v) %*% Gmat %*% v))
  
}

d <- 20
k <- 8
Gmat <- compute_weights(d, k, mc.reps = 1e5)$wmat
floor(Gmat * 1e4)

Gmat[d, d]
itersolve(Gmat)
res <- gradsolve(Gmat, nits = 1000, eps = 10)
res$v
max(res$obf)

