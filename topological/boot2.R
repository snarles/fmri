####
##  Bootstrapping testing of RSA// procrustes regression
####

library(parallel)
MCC <- 3
source("topological/rsa_boot_source.R")

p_value_dist <- function(A, B, SigmaX, SigmaY, nX, nY, theta, boot.reps, outer.reps = 100) {
  sampler <- regression_data_model_(A, B, SigmaX, SigmaY)
  pvs <- mclapply(1:outer.reps, 
                function(i) {
                  res <- sampler(nX, nY)
                  inverse_bca_test(res, theta, boot.reps)
                }, mc.cores = MCC)
  unlist(pvs)
}

pplot <- function(pvs) {
  plot(1:length(pvs)/length(pvs), sort(pvs), ylim = c(0, 1), type = "l"); abline(0, 1)  
}

p <- 5
q <- 2
sigma <- 20
nX <- 30
nY <- 30
A_0 <- randn(p, q); B_0 <- svd(randn(p, p))$u %*% A_0

f2(A_0)

pvs <- p_value_dist(A_0, B_0, eye(q)/sigma, eye(q)/sigma, nX, nY, stat.T, 1e3, 100)
pplot(pvs)
