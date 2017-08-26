

library(pracma)
source("idloss/mi_source.R")

mi_to_cor <- function(mi) {
  1/sqrt(1 + (exp(2 * mi) - 1)^(-1))
}

p <- 1
ns <- (2:200)^2
sgma <- 2# or 4.8
Xall <- randn(max(ns), p)
Yall <- Xall + sgma * randn(max(ns), p)

plot(Xall[1:1000, ]/sd(Xall), Yall[1:1000, ]/sd(Yall), asp = 1, axes = FALSE, ann = FALSE); abline(v=0); abline(h = 0)

(cor_true <- 1/sqrt(1 + sgma^2))
(mi_true <- p/2 * log(1 + sgma^(-2)))

1/sqrt(1 + (exp(2 * mi_true) - 1)^(-1))

cor_ests <- numeric()
mi_ests <- numeric()

for (ind in 1:length(ns)) {
  n <- ns[ind]
  X <- Xall[1:n, , drop = FALSE]
  Y <- Yall[1:n, , drop = FALSE]
  (cor_ests[ind] <- cor(X, Y)[1])
  (mi_ests[ind] <- nn_mi(X, Y, k=1)[1])
}

plot(sqrt(ns), cor_ests, type = "l", ylim = c(0, 1), col = "white",
     axes = FALSE, xlab = "N", ylab = "estimate")
axis(2, at = 0:5/5)
axis(1, at = 50 * 0:4, labels = (50 * 0:4)^2)
abline(cor_true, 0, lty = 2, lwd = 2, col = "grey")
lines(sqrt(ns), mi_to_cor(mi_ests), col = "red")
lines(sqrt(ns), cor_ests, col = "blue")
legend(150, 0.9, legend = c("cor", "info"), col = c("blue", "red"), lwd = 1)
