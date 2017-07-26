source("approximation/gaussian_identity_finsam.R")

p <- 10
sigma2s <- seq(0.3, 0.7, 0.1)^2
usS <- list()

for (si in 1:length(sigma2s)) {
  sigma2 <- sigma2s[si]
  sigma2_tr <- sigma2

  us <- sample_u_fs(p, sigma2, sigma2_tr, 1e4, log.p = TRUE)
  usS[[si]] <- us
}


#   K <- 2000
#   mc.reps <- 2
#   
#   t1 <- proc.time()
#   
#   mcs <- mclapply(1:mc.reps,
#                   function(i) {
#                     mus <- randn(K, p)
#                     ys <- mus + sqrt(sigma2) * randn(K, p)
#                     mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
#                     pmat <- -pdist2(mu_hats, ys)
#                     get_sub_errs(pmat, 1:K, 1:K)
#                   },
#                   mc.cores = 2)
#   
#   proc.time() - t1
#   
#   mcs2 <- do.call(cbind, mcs)
#   
#   accs <- 1-rowMeans(mcs2)
#  # plot(accs, ylim = c(0, 1))
#   
#   us <- sample_u_fs(p, sigma2, sigma2_tr, 1e4, log.p = TRUE)
#   accs2 <- sapply(1:K, function(k) mean(exp((k-1)*us)))
#   
# #  lines(accs2, col = "red")
#   
# }

#hist(exp(us))
#hist(us)
#plot(ecdf(us))
pdf("approximation/fig_mgs2.pdf")
plot(ecdf(exp(usS[[1]])), xlim = c(0,1), main = "", xlab = "u", ylab = "D(u)", col = 1)
for (i in 2:length(sigma2s)) {
  plot(ecdf(exp(usS[[i]])), col = i, lty = i, add = TRUE)
}
legend(0, 0.8, col = 1:5, lwd = 2, legend = paste("sigma =", sqrt(sigma2s)))
dev.off()


