source('fingerprint/mi_est_pipeline.R')




######
## Computation of mutual information
######

xlogx <- function(x) {
  ifelse(x < 1e-30, yes = 0, no = x * log(x))
}

# xlogx <- function(x, n.terms = 20000, threshold = 0.001) {
#   ss <- 0 * x[x < threshold]
#   for (k in 1:n.terms) {
#     ss <- ss + ((-1)^k * (-1 + x[x < threshold])^k)/k
#   }
#   y <- 0 * x
#   y[x < threshold] <- -x[x < threshold] * ss
#   y[x > threshold] = x[x > threshold] * log(x[x > threshold])
#   return(y)
# }
#plot(1:1000/100000, xlogx(1:1000/100000, n.terms = 20000, threshold = 0.001), pch = '.')


K <- 20 # number of clusters
p <- 5 # dimensionality
mus <- 40 * randn(K, p)
## Increase mc.reps until the answer below is stable
mc.reps <- 10000
y_mc <- mus[sample(K, mc.reps, TRUE), ] + randn(mc.reps, p)
dens <- exp(-pdist2(y_mc, mus)^2/2)
dens <- dens/rowSums(dens)
hx_y <- rowSums(-xlogx(dens))
(mean(hx_y))
(mi_true <- log(K) - mean(hx_y))


#####
## Estimation of MI through fingerprinting
#####


n <- 20 # number of obs for fingerprinting
y0 <- mus[sample(K, n, TRUE), ]
y1 <- y0 + randn(n, p)
y2 <- y0 + randn(n, p)
pmat <- -pdist2(y1, y2)
mi_est_pipeline(pmat)

mi_true


####
##  Do this for multiple values
####

run_gmm_experiment <- function(K, p, expand, n, mc.reps = 10000) {
  mus <- expand * randn(K, p)
  y_mc <- mus[sample(K, mc.reps, TRUE), ] + randn(mc.reps, p)
  dens <- exp(-pdist2(y_mc, mus)^2/2)
  dens <- dens/rowSums(dens)
  hx_y <- rowSums(-xlogx(dens))
  (mi_true <- log(K) - mean(hx_y))
  y0 <- mus[sample(K, n, TRUE), ]
  y1 <- y0 + randn(n, p)
  y2 <- y0 + randn(n, p)
  pmat <- -pdist2(y1, y2)
  result <- mi_est_pipeline(pmat)
  result[['mi_true']] <- mi_true
  result
}

run_simplex_experiment <- function(K, expand, n, mc.reps = 10000) {
  mus <- expand * eye(K)
  p <- K
  y_mc <- mus[sample(K, mc.reps, TRUE), ] + randn(mc.reps, p)
  dens <- exp(-pdist2(y_mc, mus)^2/2)
  dens <- dens/rowSums(dens)
  hx_y <- rowSums(-xlogx(dens))
  (mi_true <- log(K) - mean(hx_y))
  y0 <- mus[sample(K, n, TRUE), ]
  y1 <- y0 + randn(n, p)
  y2 <- y0 + randn(n, p)
  pmat <- -pdist2(y1, y2)
  result <- mi_est_pipeline2(pmat)
  result[['mi_true']] <- mi_true
  result
}

run_hcube_experiment <- function(base, order, expand, n, mc.reps = 10000) {
  K <- base^order
  arr <- array(1, dim = rep(base, order))
  mus <- expand * which(arr==1, arr.ind=TRUE) - 1
  p <- order
  y_mc <- mus[sample(K, mc.reps, TRUE), ] + randn(mc.reps, p)
  dens <- exp(-pdist2(y_mc, mus)^2/2)
  dens <- dens/rowSums(dens)
  hx_y <- rowSums(-xlogx(dens))
  (mi_true <- log(K) - mean(hx_y))
  y0 <- mus[sample(K, n, TRUE), ]
  y1 <- y0 + randn(n, p)
  y2 <- y0 + randn(n, p)
  pmat <- -pdist2(y1, y2)
  result <- mi_est_pipeline2(pmat)
  result[['mi_true']] <- mi_true
  result
}


# run_gmm_experiment(K=20, p=3, expand=20, n=20)


# run_simplex_experiment(K=200, expand=15, n=30)

run_hcube_experiment(base=6, order=4, expand=50, n=150)
