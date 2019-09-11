library(pracma)
library(lineId)
source('extrapolation/kay_method2.R')





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


K <- 10 # number of clusters
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
mi_est_pipeline <- function(pmat) {
  K <- dim(pmat)[1]
  empirical_acc <- 1- resample_misclassification(pmat, 1:K)
  
  
  acc_crit <- 0.75
  # extrapolate until 0.75 acc
  kgrid <- floor(exp((1:300)/5))
  accs_kay <- kernel_extrap3(pmat, kgrid)
  #accs_kay
  #plot(accs_kay)
  
  kchosen <- min(kgrid[accs_kay < acc_crit])
  acc_at_k <- max(accs_kay[accs_kay < acc_crit])
  
  ## get MI estimate
  #mi_est <- Ihat_LI(1-acc_at_k, kchosen)
  mi_grid <- log(kchosen) * (1:100)/200
  pi_ests <- piK(mi_grid, kchosen)
  mi_est <- max(mi_grid[(1-pi_ests) < acc_crit ])
  
  list(empirical_acc=empirical_acc, kchosen=kchosen, acc_at_k=acc_at_k, 
       mi_est = mi_est, acc_at_est = 1-piK(mi_est, kchosen))
}

n <- 10 # number of obs for fingerprinting
y0 <- mus[sample(K, n, TRUE), ]
y1 <- y0 + randn(n, p)
y2 <- y0 + randn(n, p)
pmat <- -pdist2(y1, y2)
mi_est_pipeline(pmat)

