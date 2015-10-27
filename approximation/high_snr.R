source("approximation//high_snr_source.R")
library(pracma)
library(magrittr)

p <- 10; L <- 20
pp <- seq(0.01, .99, by = 0.01)
Sigma <- cov(randn(20, p))
y <- rnorm(p)
reso <- 100
qs <- true_q(Sigma, y, pp, L)
plot(qs, guess_q(Sigma,y,pp,L))

qs[50]
pp <- pp[50]
