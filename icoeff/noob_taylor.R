## get polynomial (not Taylor) approximations a "noob" way

f <- function(x) x*log(x)

deg <- 12
xs <- seq(0, 1, 1e-4)
basis <- xs %*% t(rep(1, deg))
basis <- (basis) ^ col(basis)
ys <- f(xs); ys[1] <- 0

res <- lm(ys ~ basis)
plot(xs, ys, type = "l"); lines(xs, res$fitted.values, col = "red")
sum(is.na(res$coefficients))

#saveRDS(res$coefficients, file = "icoeff/coeffs_xlogx.rds")
