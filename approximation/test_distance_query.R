library(parallel)
library(pracma)
mcc <- 4


## query points x_i which are within distance r of x

p <- 10
n <- 1e5
xs <- randn(n, p)
nk <- floor(sqrt(n))

t1 <- proc.time()
km <- kmeans(xs, nk, iter.max = floor(nk/2))
(km_time <- proc.time() - t1)

c(max(km$size), n/nk)

## sort data by cluster (optional)
# cl_z <- km$cluster * 2 + runif(n)
# xs <- xs[order(cl_z), ]

## get maximum distances within each cluster
t1 <- proc.time()
clust_radii <- unlist(mclapply(1:nk,
                  function(clind) {
                    xsub <- xs[km$cluster == clind, ]
                    cnt <- km$centers[clind, ]
                    max(sqrt(colSums((t(xsub) - cnt)^2)))
                  },
                  mc.cores = mcc))
(dist_time <- proc.time() - t1)

hist(clust_radii)

dquery <- function(x, r, naive = FALSE) {
  if (naive) {
    alldist <- sqrt(colSums((t(xs) - x)^2))
    return(which(alldist < r))
  }
  cdist <- sqrt(colSums((t(km$centers) - x)^2))
  cdist_min <- cdist - clust_radii
  (chits <- which(cdist_min < r))
}

