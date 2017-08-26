####
## Illustration of ML vs EB
####

library(pracma)
library(magrittr)
source('transfer/source.R')
set.seed(0)
X0 <- cbind(1, randn(3, 5))
Y0 <- rbind(c(0.5, 1), c(-1, 0.1), c(.5, -1))
B <- solve(t(X0) %*% X0 + 0.001 * diag(rep(1, dim(X0)[2])), t(X0) %*% Y0)
set.seed(8)
X <- cbind(1, randn(3, 5))
Y <- X %*% B
p <- dim(Y)[2]
Sigma_e <- rbind(c(1, -.7), c(-.7, 2))
Sigma_B <- 1000 * diag(rep(1, length(B)))
Sigma_vec_B <- 10 * solve(solve(Sigma_e) %x% (t(X) %*% X) + diag(1/diag(Sigma_B)))
covs <- list(3)
for (i in 1:3) {
  covs[[i]] <- (diag(rep(1, p)) %x% t(X[i, ])) %*% 
    Sigma_vec_B %*% (diag(rep(1, p)) %x% t(t(X[i, ]))) + Sigma_e
}
svd(Sigma_e)$d
ISigma_e <- solve(Sigma_e)
ASigma_e <- isqrtm(Sigma_e)
LSigma_e <- sqrtm(Sigma_e)
tz <- 1:100/100
circ <- cbind(cos(tz * 2 * pi), sin(tz * 2 * pi))
circ2 <- circ %*% ASigma_e
temp <- list()
for (i in 1:3) temp[[paste0("mu", i)]] <- Y[i, ]
for (i in 1:2) { for (j in (i+1):3) {
  yI <- Y[i, ]
  yJ <- Y[j, ]
  temp[[paste0("m", i, j)]] <- (yI + yJ)/2
  v0 <- ((yI - yJ) %*% ISigma_e) %>% {./Norm(.)}
  vIJ <- c(v0[2], -v0[1])
  temp[[paste0("v", i, j)]] <- vIJ
}}
for (i in 1:3) {
  j <- setdiff(1:3, i)[1]
  k <- setdiff(1:3, i)[2]
  v1 <- temp[[paste0("v", min(i, j), max(i, j))]]
  v2 <- temp[[paste0("v", min(i, k), max(i, k))]]
  m1 <- temp[[paste0("m", min(i, j), max(i, j))]]
  m2 <- temp[[paste0("m", min(i, k), max(i, k))]]
  center <- m2 + (solve(cbind(v1, v2)) %*% (m1 - m2))[2]*v2  
  temp[[paste0("c", i)]] <- center
}

attach(temp)
n <- 3

pdf("poster/illus1_A.pdf")
plot(Y, pch = ".", axes = FALSE, xlab = "", ylab = "",
     xlim = range(Y[, 1]) + c(-.3, .3), ylim = range(Y[, 2]) + c(-.3, .5))
points(Y, cex = 3)
for (i in 1:3) points(Y[i, , drop = FALSE], pch = paste(i))
for (i in 1:3) lines(t(t(0.3 * circ2) + Y[i, ]), lty = 2)
for (i in 1:2) { for (j in (i+1):3) {
  x1 <- center
  x2 <- (Y[i, ] + Y[j, ])/2
  x3 <- 100 * (x2 - x1) + x1
  x4 <- -100 * (x2 - x1) + x1
  lines(rbind(x1, x3))
}}
dev.off()

detach(temp)

## have to add noise to covs
covs0 <- covs

covs <- covs0
set.seed(3)
for (i in 1:3) {
  covs[[i]] <- 0.7 * covs[[i]] + 2 * (randn(2,2) %>% {t(.)%*%.})
}
covs
## find points with equal prob of being in classes

trdist <- function(x) {
  ans <- numeric(3)
  for (i in 1:3) ans[i] <- t(x- Y[i, ]) %*% solve(covs[[i]]) %*% (x - Y[i, ])
  ans
}
discrp <- function(x) var(trdist(x))

center <- optim(colMeans(Y), discrp, method = "BFGS")$par

js <- c(2, 3, 1)

discrps <- function(x) {
  ans0 <- trdist(x)
  ans <- numeric(3)
  for (i in 1:3) ans[i] <- (ans0[i] - ans0[js[i]])^2
  ans
}

pdf("poster/illus1_B.pdf")
plot(Y, pch = ".", axes = FALSE, xlab = "", ylab = "",
     xlim = range(Y[, 1]) + c(-.3, .3), ylim = range(Y[, 2]) + c(-.3, .5))
for (i in 1:3) points(Y[i, , drop = FALSE], pch = paste(i))
lvdists <- 1:6/7
for (i in 1:3) {
  circ3 <- circ %*% isqrtm(covs[[i]])
  for (lv in lvdists) {
    lines(t(t(lv * circ3) + Y[i, ]), lty = 2)
  }
}  
sns <- c(1, -1, 1)
for (i in 1:3) {
  j <- js[i]
  m <- (Y[i, ] + Y[j, ])/2
  v <- (m - center) %>% {./Norm(.)}
  path <- t(t(rep(1, 21))) %*% center + sns[i] * t(t(0:20/10)) %*% v
  res <- apply(path, 1, function(x) {
    optim(x, function(x) {
      discrps(x)[i]
    }, method = "BFGS")$par
  })
  lines(t(res))
}
dev.off()