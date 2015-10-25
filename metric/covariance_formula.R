library(pracma)

TR <- function(A) sum(diag(A))
cov_formula <- function(Sigma1, Sigma2) 2 * TR(Sigma1 %*% Sigma2)

Sigma1 <- cov(randn(10, 5))
Sigma2 <- cov(randn(10, 5))
nms <- do.call(rbind,lapply(1:1e3, function(i) {
  z <- rnorm(5)
  c(t(z) %*% Sigma1 %*% z, t(z) %*% Sigma2 %*% z)}))

cov(nms[, 1], nms[, 2])
cov_formula(Sigma1, Sigma2)
