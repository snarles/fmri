source("approximation//high_snr_source.R")
library(pracma)
library(magrittr)

p <- 10; L <- 1000
pp <- seq(0.01, .99, by = 0.01)
Sigma <- cov(randn(1000, p)) + 20 * eye(p)
y <- 20 * rnorm(p)
reso <- 100
xs <- true_q(Sigma, y, pp, L)
xdowns <- find_x_bar(Sigma, y, pp, L)$xdown
plot(xs, xdowns)
sum(xs > xdowns)
median(xs/xdowns)
