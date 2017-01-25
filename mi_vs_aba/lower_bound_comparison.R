## compare lower bound to specific results

source("approximation/gaussian_identity.R")

library(lineId)


K <- 5

Is <- seq(0, 10, by = 0.1)
risks_inf <- lineId::piK(sqrt(2 * Is), K)

risks_lb <- seq(0, 1-(1/K), by = 0.01)[-1]
Is_lb <- sapply(risks_lb, function(x) aba_to_mi_lower(K, 1 - x))



sigs <- exp(seq(-8.5, 4, by = 0.5))
risks_2 <- sapply(sigs, function(x) mc_ident2(1, x, K, 3e5))
rhos <- 1/sqrt(1 + sigs)
Is_2 <- -1/2 * log(1 - rhos^2)

par(bg = "grey")

plot(NA, NA, xlab = "IdRisk", ylab = "MI", ylim = c(0, 4), xlim = c(0, 1))
lines(risks_lb, Is_lb, lwd = 2, col = "black")
lines(risks_inf, Is, lwd = 2, col = "blue")
lines(risks_2, Is_2, lwd = 2, col = "yellow")
legend(0.7, 4, c("bound", "p=Inf", "p=2"), col = c("black", "blue", "yellow"),
       lwd = 2)
title("k = 5")
