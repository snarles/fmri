## Lambda(beta) = int exp(beta * t^(k-1)) dt
### claim: (k-1) Lambda(beta) converges as k to infty

delt <- 1e-4
ts <- seq(0, 1, by = delt)
ks <- 1:200
bts <- exp(seq(-2, 3, by = 0.01))
#curves <- sapply(ks, function(k) delt * rowSums(exp(bts %*% t(ts)^(k-1))))
curves <- sapply(ks, function(k) {
  (k-1) * (delt * rowSums(exp(bts %*% t(ts)^(k-1))) - 1)
})

is <- 2:100
limit_curve <- sapply(bts, function(bt) {
  sum(exp(is * log(bt) - lgamma(is + 1) - log(is)))
})
matplot(bts, curves[, 5 * 1:10], type = "l")

plot(ks, curves[4, ], type = "l")

plot(bts, curves[, 200], type = "l")
lines(bts, limit_curve, col = "red")

