### fit exponentially weighted mixture

library(nnls)


expmix <- function(ws, as, xs) {
  as.numeric(ws %*% exp(t(t(as)) %*%  t(xs)))
}

expbasis <- function(as, xs) {
  t(exp(t(t(as)) %*%  t(xs)))
}

fit_expmix <- function(as, y) {
  X <- expbasis(as, xs)
  res <- nnls(X, y)
  sol <- res$x
  sol[sol < 1e-10] <- 0
  fit_a <- as[sol > 0]
  fit_w <- sol[sol > 0]
  list(a = fit_a, w = fit_w)
}

# as <- -(1:3)
# xs <- seq(0, 10, 0.1)
# ws <- rep(1, 3)/3
# basis_a <- -seq(0, 4, 0.1)
# y <- expmix(ws, as, xs)
# X <- expbasis(basis_a, xs)
# res <- nnls(X, y)
# res$x
# fit_expmix(basis_a, y)
