####
##  Bootstrapping testing of RSA// procrustes regression
####

library(parallel)
MCC <- 39
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

pplot <- function(pvs, ...) {
  plot(1:length(pvs)/length(pvs), sort(pvs), ylim = c(0, 1), type = "l", ...); abline(0, 1)  
}

p <- 20 # number of neurons
q <- 2 # dimension of regressor
nX <- 50
nY <- 50
boot.reps <- 1e3
o.reps <- 5e2

f2a <- 10
f2d <- 5
A_0 <- randn(p, q)
DD <- randn(p, q); DD <- sqrt(f2d) * DD/sqrt(f2(DD))
A_0 <- sqrt(f2a) * A_0/sqrt(f2(A_0))
B_0 <- svd(randn(p, p))$u %*% A_0
B_1 <- B_0 + DD

## Null case

pvs <- p_value_dist(A_0, B_0, eye(q), eye(q), nX, nY, stat.T, boot.reps, o.reps)
pplot(pvs, ann = FALSE)
title(paste("h0 T, p=",p,"q=",q,"nX=",nX,"nY=",nY,"B=",boot.reps),
      sub = paste("f2a =", f2a))

pvs <- p_value_dist(A_0, B_0, eye(q), eye(q), nX, nY, stat.S, boot.reps, o.reps)
pplot(pvs, ann = FALSE)
title(paste("h0 S, p=",p,"q=",q,"nX=",nX,"nY=",nY,"B=",boot.reps),
      sub = paste("f2a =", f2a))

## Alt case

pvs <- p_value_dist(A_0, B_1, eye(q), eye(q), nX, nY, stat.T, boot.reps, o.reps)
pplot(pvs, ann = FALSE)
title(paste("h1 T, p=",p,"q=",q,"nX=",nX,"nY=",nY,"B=",boot.reps),
      sub = paste("f2a =", f2a, "f2d=", f2d))

pvs <- p_value_dist(A_0, B_1, eye(q), eye(q), nX, nY, stat.S, boot.reps, o.reps)
pplot(pvs, ann = FALSE)
title(paste("h1 S, p=",p,"q=",q,"nX=",nX,"nY=",nY,"B=",boot.reps),
      sub = paste("f2a =", f2a, "f2d=", f2d))

