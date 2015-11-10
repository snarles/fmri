library(pracma)
library(lineId)
library(magrittr)

## minimize ||HG - J||^2 subject to G orthogonal
rotmin <- function(H, J) {
  HtJ <- t(H) %*% J
  isqrtm(HtJ %*% t(HtJ)) %*% HtJ
}


### TESTS ###
# p <- 10
# H <- randn(p)
# J <- randn(p)
# G <- rotmin(H, J)
# f2(H%*%G, J)
# min(sapply(1:1000, function(i) f2(H%*%svd(randn(p))$u, J)))
# f2(G %*%t(G), eye(p))
