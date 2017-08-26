#### Do simulations k = 2 to 50
source("extrapolation/simulation_source.R")
K <- 50
synth_data <- gen.data(p = 10, sigma = 1, K = K, r1 = 20, r2 = 30)
tabs <- list()
ks <- 3:50
for (k in ks) {
  tab <- build_extrapolation_table(synth_data, k, K)
  tabs <- c(tabs, list(tab))
}

saveRDS(tabs, file = "extrapolation/simres.rds")

colnames(tab)

cols <- c("black", "green", "red", "blue")
WW <- 4
HH <- 4

pdf("extrapolation/sim_qda.pdf", width = WW, height = HH)
method <- "qda"
accs <- sapply(tabs, function(v) v[method, c(1, 3, 7, 8)])
matplot(ks, t(accs), type = "l", lwd = 3, col = cols, lty = 1, ylim = c(0, 1))
abline(accs["acc_sub", length(ks)], 0, lty = 2, lwd = 2)
title("QDA")
dev.off()

pdf("extrapolation/sim_glmnet.pdf", width = WW, height = HH)
method <- "glmnet"
accs <- sapply(tabs, function(v) v[method, c(1, 3, 7, 8)])
matplot(ks, t(accs), type = "l", lwd = 3, col = cols, lty = 1, ylim = c(0, 1))
abline(accs["acc_sub", length(ks)], 0, lty = 2, lwd = 2)
title("Logistic")
dev.off()

pdf("extrapolation/sim_nnet.pdf", width = WW, height = HH)
method <- "nnet"
accs <- sapply(tabs, function(v) v[method, c(1, 3, 7, 8)])
matplot(ks, t(accs), type = "l", lwd = 3, col = cols, lty = 1, ylim = c(0, 1))
abline(accs["acc_sub", length(ks)], 0, lty = 2, lwd = 2)
title("Neural Net")
dev.off()
