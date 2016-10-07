## Dependence, distance, and information

library(MASS)

rho <- 0
rho <- 0.5
col1 <- rgb(1, 0, 0, 1)
col2 <- rgb(0, 0, 1, 1)

plot1 <- function(rho) {
  layout(1)
  set.seed(0)
  mat <- matrix(c(1, rho, rho, 1), 2, 2)
  xy <- mvrnorm(2000, mu = c(0, 0), Sigma = mat)
  plot(NA, NA, xlim = c(-4, 4), ylim = c(-3, 3), xlab = "x", ylab = "y")
  points(xy[xy[, 1] >= 0, ], col = col1, cex = 0.3)
  points(xy[xy[, 1] < 0, ], col = col2, cex = 0.3)
}

plot2 <- function(rho) {
  #layout(matrix(1:2), 2, 1)
  set.seed(0)
  mat <- matrix(c(1, rho, rho, 1), 2, 2)
  xy <- mvrnorm(2000, mu = c(0, 0), Sigma = mat)
  x <- xy[, 1]
  y <- xy[, 2]
  bks <- seq(-4, 4, by = 0.2)
  hist(y[x >= 0], breaks = bks, xlim = c(-4, 4), col = col1, main = "Y|X>0", xlab = "x", ylim = c(0, 150))
  #hist(x[x <= 0 & y >= 0], breaks = bks, xlim = c(-4, 4), col = col2, add = TRUE)
  hist(y[x <= 0], breaks = bks, xlim = c(-4, 4), col = col2, main = "Y|X<0", xlab = "x", ylim = c(0, 150))
  #hist(x[x <= 0 & y <= 0], breaks = bks, xlim = c(-4, 4), col = col2, add = TRUE)
  print(list(acc = mean(x * y >= 0)))
}

rho <- 0.9

classific <- function(rho) {
  mat <- matrix(c(1, rho, rho, 1), 2, 2)
  xy <- mvrnorm(1e5, mu = c(0, 0), Sigma = mat)
  x <- xy[, 1]
  y <- xy[, 2]
  z <- (x >= 0)
  res <- glm(z ~ poly(y, degree = 3), family = "binomial")
  plot(y[order(y)], res$fitted.values[order(y)], type = "o")
  fv <- res$fitted.values
  mean(pmax(fv, 1-fv))
}

library(imputeTS)

kldiv <- function(rho, reso = 1e6) {
  mat <- matrix(c(1, rho, rho, 1), 2, 2)
  xy <- mvrnorm(reso, mu = c(0, 0), Sigma = mat)
  x <- xy[, 1]
  y <- xy[, 2]
  d1 <- density(x[y > 0], from = -5, to = 5, n = 1000)
  d2 <- density(x[y < 0], from = -5, to = 5, n = 1000)
  delt <- d1$x[2] - d1$x[1]
  ss <- (d1$y * log(d1$y/d2$y))
  filt <- is.nan(ss) | is.infinite(ss)
  #plot(df$x, res$fitted, type = "l")
  if (sum(filt) > 0) {
    ss[filt] <- NA
    df <- data.frame(y = ss[!filt], x = d1$x[!filt])
    res <- loess(y ~ x, data = df, span = 0.1)
    mvs <- predict(res, data.frame(x = d1$x[filt]))
    ss[filt] <- mvs
  }
  #plot(d1$x, ss, type = "l")
  # points(d1$x[is.na(ss)], rep(0, sum(is.na(ss))) , col = "red")
  sum(ss, na.rm = TRUE) * delt
}

plot1(0)
plot1(0.8)

plot2(0)
plot2(0.8)

kldiv(0)
kldiv(0.1)
kldiv(0.2)
kldiv(0.3)
kldiv(0.8, 1e8)

rhos <- seq(0, 0.95, by = 0.05)
kls <- sapply(rhos, kldiv)
layout(1)
plot(rhos, kls, type = "l")

classific(0)
classific(0.1)
classific(0.2)
classific(0.3)
classific(0.4)
classific(0.8)
