## Taylor expansion question
## what is the optimal u0
## minimizing int_0^1 (u-u0)^d u^{k-2} du

delta <- 0.01
us <- seq(0, 1, delta)
d <- 4 # needs to be even
k <- 50
plot(us, us^(k-2), type = "l")
ints <- sapply(us, function(u0) delta * sum((us-u0)^d * us^(k-2)))

plot(us, log(ints), type = "l")
coefs <- choose(d, 0:d) * (-1)^(0:d)/(k - 1 + rev(0:d))
lines(us, log(sapply(us, function(u) sum(u^(0:d) * coefs))), col = "red")
abline(v = (k-1)/k)
