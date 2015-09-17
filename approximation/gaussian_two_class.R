####
##  Error for random two-class gaussian classification
####



source("approximation//gaussian_2c_source.R")

Sigma <- 0.05 * rchisq(1, 5) * cov(randn(40, 10))
mc2c(Sigma, 1e4)
mc2c_1(Sigma, 1e4)
mc2c_2(Sigma, 1e4)
mc2c_3(Sigma, 1e4)

p <- length(as)
lambdas <- eigen(Sigma)$values
as <- (2 * lambdas)^2
zs <- randn(p, mc.reps)
ss <- colSums(as * (zs^2))
c(mean(ss), sum(as))
c(var(ss), 2 * sum(as^2))




build_mc2c_table(eye(2), 1/4, returnAll = TRUE)

### Calibration record
## 
## build_mc2c_table(eye(2), 1/4) == 0.6660156
## build_mc2c_table(eye(3), 1/4) == 0.3896484
## build_mc2c_table(eye(4), 1/4) == 0.2744141
## build_mc2c_table(eye(5), 1/4) == 0.2114258



scalings <- unlist(mclapply(1:20,
                            function(k) build_mc2c_table(eye(k), 1/4), mc.cores = 3))

v2copypasta(scalings)


matA <- diag(c(1, rep(0, 9)))
matB <- eye(10)


fracs <- c(0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, .5)
fracs <- c(fracs, 1 - rev(fracs)[-1])
Sigmas <- lapply(fracs, function(x) x * matA + (1-x) * matB)
scalings <- unlist(mclapply(Sigmas,
                            function(s) build_mc2c_table(s, 1/4),
                            mc.cores = 7))
v2copypasta(scalings)
