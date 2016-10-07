library(pracma)
p <- 6
n <- 1000

X <- randn(n, p)
z <- (rowSums(X) > 0) + 0

ngrid <- qnorm(seq(0.001, 0.999, by = 0.001))
mis <- numeric()
cas <- numeric()
for (i in 1:5) {
  vfrac <- sqrt(i/(p-i))
  ps <- pnorm(vfrac * ngrid)
  ce <- -mean(ps * log(ps) + (1-ps) * log(1-ps))
  mis[i] <- log(2) - ce
  cas[i] <- mean(pmax(ps, 1-ps))
}

plot(mis, type = "o", xlab = "k", ylab = "", lwd = 2,
     main = "Mutual information", ylim = c(0, 0.45))
plot(cas, type = "o", xlab = "k", ylab = "", lwd = 2,
     main = "Bayes accuracy", ylim = c(0.5, 1))
