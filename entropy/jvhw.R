library(R.matlab)
library(pracma)

fingerprint <- function(freqs) {
  tt <- table(freqs)
  vs <- as.numeric(names(tt))
  ff <- numeric(max(vs))
  ff[vs] <- tt
  ff
}

h_mle <- function(freqs) {
  ps <- freqs/sum(freqs)
  -sum(ps * log(ps))
}

Vconst <- c(0.3303, -0.3295, 0.4679)
poly_coeff_r <- readMat('entropy//poly_coeff_r.mat')

h_jvhw <- function(freqs) {
  n <- sum(freqs)
  order <- 4 + ceiling(1.2*log(n))
  coeff <- poly_coeff_r[[1]][order][[1]][[1]]
  ff <- fingerprint(freqs)
  prob <- 1:length(ff)/n
  c_1 <- 0
  if (n > 15 && ff[1] > 0) {
    c_1 <- sum(Vconst * c(log(n), log(ff[1]), 1))
    c_1 <- pmax(c_1, 1/(1.9*log(n)))
  }
  c_1
  x <- prob
  g_coeff <- coeff
  K <- length(g_coeff) - 1
  thres = 4*c_1*log(n)/n
  output = zeros(length(x), length(c_1))
  region_large = x>thres
  region_nonlarge = !region_large
  region_mid = x>thres/2 & region_nonlarge
  output[region_large] = -x[region_large] * log(x[region_large]) + 1/(2*n); 
  if (sum(region_nonlarge) > 0) {
    x1 <- x[region_nonlarge]
    q = 0:(K-1)
    n * x1
    temp0 <- cbind(thres, (repmat(n * t(t(x1)), 1, K) - repmat(q, length(x1), 1))/(thres * repmat(n - q, length(x1), 1)))
    temp <- t(apply(temp0, 1, cumprod))
    output[region_nonlarge] <- temp %*% t(g_coeff) - x1 * log(thres)    
  }
  ratio = 2*x[region_mid]/thres - 1
  output[region_mid] = ratio*(-x[region_mid] * log(x[region_mid]) + 1/(2*n)) + 
    (1-ratio)*output[region_mid]
  est = sum(ff*output)
  est
}
