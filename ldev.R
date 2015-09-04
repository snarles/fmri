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

####
##  Testing Bayesian methods
####

#' Tests Bayesian method
#'
#' Let (B, Y) have a joint distribution, where Y is observed
#' Then for any functions f(B), g(Y) we have
#' E[f(B)g(Y)] = E[h(Y)g(Y)]
#' where h(Y) = E[f(B)|Y]

ntrial <- 1000

#' Normal cf
charfn <- function(tt, mu, Sigma) exp(1i * sum(mu * tt) - (t(tt) %*% Sigma %*% tt)/2)

n <- 100; pX <- 20; pY <- 30
X <- randn(n, pX)
Sigma_b <- diag(rexp(pY))
Sigma_e <- cor(randn(3 * pY, pY))
Sigma_t <- cor(randn(3 * n, n))
sst <- sqrtm(Sigma_t)
sse <- sqrtm(Sigma_e)
Sigma_B <- (eye(pX) %x% Sigma_b)
Sigma_Y <- (X %*% t(X)) %x% Sigma_b + Sigma_t %x% Sigma_e
Sigma_BY <- rbind(cbind(Sigma_B, t(X) %x% Sigma_b),
                  cbind(X %x% Sigma_b, Sigma_Y))

BYs <- lclapply(1:ntrial, function(i) {
    B <- t(sqrt(diag(Sigma_b)) * randn(pY, pX))
    E <- sst %*% randn(n, pY) %*% sse
    list(B = B, Y = X %*% B + E)
  }, mc.cores = 3)

fmu <- rnorm(pX * pY)
nfmu <- (t(fmu) %*% Sigma_B %*% fmu)[1]
fmu <- fmu/sqrt(nfmu) * runif(1)

gmu <- rnorm(n * pY)
ngmu <- (t(gmu) %*% Sigma_Y %*% gmu)[1]
gmu <- gmu/sqrt(ngmu) * runif(1)

fBs <- sapply(BYs$B, function(B) exp(1i * sum(as.numeric(B) * fmu)))
gYs <- sapply(BYs$Y, function(Y) exp(1i * sum(as.numeric(Y) * gmu)))
fBgYs <- fBs * gYs

fBs0 <- charfn(fmu, rep(0, pX * pY), Sigma_B)[1]
gYs0 <- charfn(gmu, rep(0, n * pY), Sigma_Y)[1]
fBgYs0 <- charfn(c(fmu, gmu), rep(0, (n + pX) * pY), Sigma_BY)[1]

## c(empirical value, theoretical value)
c(mean(fBs), fBs0)
c(mean(gYs), gYs0)
c(mean(fBgYs), fBgYs0)

## computation using h
Y <- BYs$Y[[1]]
hY <- function(Y) {
  res <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, TRUE)
  charfn(fmu, as.numeric(res$Mu), res$Cov)[1]
}
hYs <- unlist(mclapply(BYs$Y, hY, mc.cores = 3))
c(mean(fBs), mean(hYs), fBs0)
c(mean(fBgYs), mean(hYs * gYs), fBgYs0)




