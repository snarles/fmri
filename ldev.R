## DEVELOPMENT OF lineId package
## temporary file to be moved later

library(devtools)
library(roxygen2)
roxygenize('lineId')
build('lineId')
install('lineId')
library(lineId)

help(tkron_d_kron)


a1 <- 5
a2 <- 3
b1 <- 3
b2 <- 4

f2
A <- pracma::randn(a1, a2)
B <- pracma::randn(b1, b2)
cc <- rnorm(a1 * a2)
ans <- tkron_d_kron(A, B, cc)
ansn <- t(A %x% B) %*% diag(cc) %*% (A %x% B)
f2(ans, ansn)

