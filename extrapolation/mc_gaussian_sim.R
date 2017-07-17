source("approximation/gaussian_identity_finsam.R")

p <- 10
sigma2 <- 1
sigma2_tr <- 0.1
K <- 2000
mc.reps <- 2

t1 <- proc.time()

mcs <- mclapply(1:mc.reps,
              function(i) {
                mus <- randn(K, p)
                ys <- mus + sqrt(sigma2) * randn(K, p)
                mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
                pmat <- -pdist2(mu_hats, ys)
                get_sub_errs(pmat, 1:K, 1:K)
              },
              mc.cores = 2)

proc.time() - t1

mcs2 <- do.call(cbind, mcs)

accs <- 1-rowMeans(mcs2)
plot(accs, ylim = c(0, 1))

us <- sample_u_fs(p, sigma2, sigma2_tr, 1e4, log.p = TRUE)
accs2 <- sapply(1:K, function(k) mean(exp((k-1)*us)))

lines(accs2, col = "red")

hist(exp(us))
hist(us)
plot(ecdf(us))
plot(ecdf(exp(us)), xlim = c(0,1), main = "")

