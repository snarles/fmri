## Extrapolation method based on Kernel Density estimation


## quickly generate some data
p <- 10
K <- 1e4
sigma2 <- 0.3
sigma2_tr <- sigma2
mus <- randn(K, p)
ys <- mus + sqrt(sigma2) * randn(K, p)
mu_hats <- mus + sqrt(sigma2_tr) * randn(K, p)
pmat <- -pdist2(ys, mu_hats)

TOL <- 1e-9
raccs <- sapply(1:K, function(ind) {
  dens <- density(pmat[ind, ], bw = "bcv", from = min(pmat[ind, ]), to = max(pmat[ind, ]) + 7 * sd(pmat[ind, ]), n = 2048)
  stopifnot(dens$y[length(dens$y)] < TOL)
  ind.at <- min(which(dens$x > pmat[ind, ind]))
  stopifnot(abs(dens$x[ind.at] - dens$y[ind.at + 1]) < TOL)
  (racc <- 1 - sum(dens$y[dens$x > pmat[ind, ind]]) * (dens$x[2] - dens$x[1]))
  racc  
})
