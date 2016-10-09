
layout(1)

k <- 5
bts <- seq(1/k, 1, by = 0.001)
iotas <- bts * log(bts) + (1-bts) * log((1-bts)/(k-1)) + log(k)
plot(iotas, bts, type = "l", ylim = c(0, 1), xlim = c(0, log(k)), lwd = 2)
max(iotas, na.rm = TRUE)
log(k)
max(bts, na.rm = TRUE)

ps <- seq(1/k, 0.999, by = 0.001)
fs <- log(k) - (1-ps)*log(k-1) + (1-ps)*log((1-ps)) + ps*log(ps)
lines(fs, ps, col = "blue", lwd = 2)
