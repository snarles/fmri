## DEVELOPMENT OF lineId package
## temporary file to be moved later

library(devtools)
library(roxygen2)
roxygenize('lineId')
load_all('lineId')
build('lineId')
install('lineId')
library(lineId)

help(tkron_d_kron)


a1 <- 5; a2 <- 3; b1 <- 3; b2 <- 4
A <- pracma::randn(a1, a2)
B <- pracma::randn(b1, b2)
cc <- rnorm(a1 * b1)
ans <- tkron_d_kron(A, B, cc)
ansn <- t(A %x% B) %*% diag(cc) %*% (A %x% B)
f2(ans, ansn)

a1 <- 5; a2 <- 3; b1 <- 3; b2 <- 4
A <- pracma::randn(a1, a2)
B <- pracma::randn(b1, b2)
cc <- rnorm(a2 * b2)
ans <- kron_v(A, B, cc)
ansn <- (A %x% B) %*% cc
f2(ans, ansn)

n <- 100; pX <- 20; pY <- 30
X <- randn(n, pX)
B <- randn(pX, pY)
Sigma_b <- eye(pY)
Sigma_e <- cor(randn(3 * pY, pY))
Sigma_t <- cor(randn(3 * n, n))
E <- sqrtm(Sigma_t) %*% randn(n, pY) %*% sqrtm(Sigma_e)
Y <- X %*% B + E
res <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, TRUE, TRUE)
f2(res$Mu, B)
f2(solve(t(X) %*% X, t(X) %*% Y), B)

