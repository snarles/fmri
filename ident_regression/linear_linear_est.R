############
## Estimation in the linear/linear case
############

## empirical misclassification rate

library(parallel)
library(class)

simulate <- function(bt, sigma2_x, sigma2_eps,
                     al, k_cl, n_trials, seed = 0, mc.cores = 25) {
  simulate0 <- function(seed) {
    set.seed(seed)
    mrs <- numeric(n_trials)
    for (i in 1:n_trials) {
      xs <- rnorm(k_cl, 0, sqrt(sigma2_x))
      ys <- bt * xs + rnorm(k_cl, 0, sqrt(sigma2_eps))
      yhats <- xs * al
      te_cl <- knn(t(t(yhats)), t(t(ys)), 1:k_cl, k = 1)
      mr <- sum(te_cl != 1:k_cl)/k_cl
      mrs[i] <- mr
    }
    return(mrs)    
  }
  if (length(seed) == 1) {
    return(simulate0(seed))
  }
  else {
    return(mclapply(seed, simulate0, mc.cores = mc.cores))
  }
}

## probability that |N(mu, sigma2)| > thres
p_folded <- function(mu, sigma2, thres, empirical = FALSE) {
  mu_conv <- mu/sqrt(sigma2)
  thres_conv <- abs(thres)/sqrt(sigma2)
  if (empirical) {
    return(sum(abs(rnorm(1e6, mu, sqrt(sigma2))) < thres)/1e6)
  }
  1 - pnorm(-mu_conv + thres_conv) + pnorm(-mu_conv - thres_conv)
}

## misclassification rate based on conditional dist
theory1 <- function(bt, sigma2_x, sigma2_eps,
                    al, k_cl, res = 30) {
  # normal distribution grid
  g0 <- -res:res
  z2 <- cbind(rep(g0, 2*res + 1),
              rep(g0, each = 2*res + 1)) /
    sqrt(res)
  d2 <- exp(-.5 * z2[,1]^2 - .5 * z2[,2]^2 )
  d2 <- d2/sum(d2)
  xs <- sqrt(sigma2_x) * z2[,1]
  eps <- sqrt(sigma2_eps) * z2[,2]
  mus <- bt * xs + eps
  sigma2 <- al^2 * sigma2_x
  thress <- (bt - al) * xs + eps
  ps <- p_folded(mus, sigma2, thress)
  1 - sum(ps^(k_cl-1) * d2)
}


bt <- 2
sigma2_x <- 1
sigma2_eps <- 1
ps <- simulate(bt, sigma2_x, sigma2_eps, 2, 3, 10000, 1:25)
ps <- unlist(ps)
mean(ps)
sd(ps)/sqrt(length(ps))
theory1(bt, sigma2_x, sigma2_eps, bt, 3, 400)

bts <- (-200:200)/20
k_cl <- 4
res <- length(bts)
iloop1 <- function(i) {
  ans <- numeric(res)
  if (bts[i] == 0) {
    ans <- rep(1 - 1/k_cl, res)
  }
  for (j in 1:res) {
    ans[j] <- theory1(bts[j], sigma2_x, sigma2_eps, bts[i], k_cl, 300)
  }
  ans
}

proc.time()
temp <- mclapply(1:res, iloop1, mc.cores = 26)
proc.time()

lines(newtemp[[4]], col = "red")

sapply(temp, length)
rmat <- do.call(cbind, temp)
rmat[, 198:204] <- do.call(cbind, newtemp)

rmat[, 201] <- 1 - (1/k_cl)

#pdf("paper/rmat.pdf")
filled.contour(bts, bts, rmat, xlab = expression(hat(beta)), ylab = expression(beta))
title(expression(paste("R(", beta, "; ", hat(beta), ")")))
#dev.off()

## Optimal beta hat given the prior
opt_bth <- function(mu, sigma2) {
  dn <- dnorm(bts, mean = mu, sd = sqrt(sigma2))
  dn <- dn/sum(dn)
  lala <- dn %*% rmat
  bth = bts[lala == min(lala)]
  bth
}

####
## SMOOTH OUT THE RMAT
####
rmat_old <- rmat
rmat <- rmat_old



nn <- dim(rmat)[1]
s2d <- diag(rep(1, nn))
s2d[abs(row(s2d) - col(s2d)) == 1] <- -.5
s2ds <- s2d %*% rmat
s2ds[1, ] <- 0; s2ds[nn, ] <- 0
plot(abs(s2ds)[, 100])
plot(abs(s2ds)[, 200])
plotit <- function(k, thres) {
  s2ds <- s2d %*% rmat
  s2ds[1, ] <- 0; s2ds[nn, ] <- 0
  lala <- which(abs(s2ds)[, k] > thres)
  plot(rmat[, k], type = "l")
  points(lala, rmat[lala, k], col = "red")
}
thres <- 0.001
smooth_rmat <- function(rmat, thres) {
  for (k in which(apply(abs(s2ds), 2, max) > thres)) {
    lala <- which(abs(s2ds)[, k] > thres)
    for (i in lala) {
      minl <- max(1, i - 10)
      maxl <- min(nn, i + 10)
      inds <- setdiff(minl:maxl, lala)
      if (length(inds) < 4) inds <- setdiff(minl:maxl, i)
      y <- rmat[inds, k]
      x <- cbind(1, inds, inds^2, inds^3)
      cf <- lm(y ~ x + 0)$ coefficients
      rmat[i, k] <- sum(c(1, i, i^2, i^3) * cf)
    }
  }
  rmat
}
rmat <- smooth_rmat(rmat, 0.01)
rmat <- smooth_rmat(rmat, 0.005)
rmat <- smooth_rmat(rmat, 0.003)
rmat <- smooth_rmat(rmat, 0.001)


plotit(210, 0.001)

plotit(203, 0.001)
plotit(202, 0.001)
plotit(201, 0.001)
lines(rmat[, 202], col = "red")
lines(rmat[, 203], col = "blue")

plot((rmat[, 200] + rmat[, 202])/2, type = "l")


plotit(200, 0.001)

plotit(190, 0.001)
plotit(150, 0.001)
plotit(100, 0.001)

sigma2 <- 1
cand_bts <- (-200:200)/40
temp1 <- function(i) {
  opt_bth(cand_bts[i], sigma2)
}
#temp1(1)
pdf("poster/zero2.pdf")
res <- sapply(1:length(cand_bts), temp1)
layout(1)
plot(cand_bts, res, type = 'l', xlab = "Point estimate", ylab = "Optimal")
lines(cand_bts, cand_bts, lty = 2)
dev.off()
#title(expression(paste(sigma^2, "=")))
#title(paste("        ", sigma2))

## minimal estimate phenomenon

min(abs(res))
## sigma2 = 10,  3.3
##          5,   2.55
##          0.5, 1.1

max(rmat)
l1n <- function(v) v/sum(v)
sds <- c(0.25, .3, .5, 1, 2, 3)
range(bts)
#plot(NA, NA, xlim = range(bts), ylim = c(0, 1),
#     xlab = expression(beta), ylab = "R")
layout(matrix(1:length(sds), length(sds), 1))
oldmar <- par("mar")
par(mar = c(0, 0, 0, 0))
for (i in 1:length(sds)) {
  y <- l1n(dnorm(bts, mean = 0, sd = sds[i])) %*% rmat
  plot(bts, y,
        type = 'l', col = rainbow(length(sds))[i])  
}
par(mar = oldmar)

plot(rmat[200, ], type = "l")
save(rmat, file = "rmat3.RData")
