## Lambda(beta) = int exp(beta * t^(k-1)) dt
### claim: (k-1) Lambda(beta) converges as k to infty

delt <- 1e-4
ts <- seq(0, 1, by = delt)
ks <- 1:200
bts <- exp(seq(-2, 3, by = 0.01))
#curves <- sapply(ks, function(k) delt * rowSums(exp(bts %*% t(ts)^(k-1))))
curves <- sapply(ks, function(k) (k-1) * delt * rowSums(exp(bts %*% t(ts)^(k-1))))

limit_curve <- colSums((t(exp(bts %*% t(ts))) * delt/ts)[-1, ])
matplot(bts, curves[, 5 * 1:10], type = "l")
lines(bts, limit_curve, type = "l", lwd = 2)
