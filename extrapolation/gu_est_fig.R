source("extrapolation/psuedolikelihood.R")

layout(1)

## Generate data
set.seed(0)
n <- 100
us <- seq(0, 1, 0.01)
gu <- cumsum(runif(length(us)) * (1:length(us))^4)
gu <- gu/sum(gu)
usamp <- sample(us, n, TRUE, prob = gu)
k <- 10
Ys <- rbinom(n, k, usamp)
(ws <- sapply(0:k, function(i) sum(Ys == i)))
(momK <- mean(binmom(Ys, k, k)))
usk <- us^k

res <- fit_pm_models(Ys, k, us)
moments_fit <- matrix(0, 50, 4)
for (i in 1:nrow(moments_fit)) {
  moments_fit[i, 1] <- sum(gu * us^i)
  moments_fit[i, 2] <- sum(res$gu_mple * us^i)
  moments_fit[i, 3] <- sum(res$gu_mono * us^i)
  moments_fit[i, 4] <- sum(res$gu_mm * us^i)
}


plot(us, gu, type = "l", ylim = c(0, 0.3), lwd = 3)
lines(us, res$gu_mple, col = "blue", lwd = 3)
#lines(us, res$gu_mono, col = "green", lwd = 3)
lines(us, res$gu_mm, col = "red", lwd = 3)
lines(us, gu, type = "l", lwd = 3)

plot(1:50, moments_fit[, 1], type = "l", lwd = 3, ylim = c(0, 1))
matplot(10:50, moments_fit[10:50, -3], type = "l", lwd = 3, col= c("black", "blue", "red"),add = TRUE,
        lty = c(1, 1, 1, 1))
abline(v = 10, lty = 2, col = "grey", lwd = 3)
