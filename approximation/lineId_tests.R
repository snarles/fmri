library(lineId)
library(pracma)

mu <- 4
sigma2 <- 8
piK(mu, 200, mc.reps = 1000, sigma2)

raw <- randn(1e5, 200)
mins <- apply(raw[, -1], 1, max)
ggg <- raw[, 1] * sqrt(sigma2) +mu
mean(ggg < mins)
