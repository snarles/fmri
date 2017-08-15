library(parallel)
library(pracma)
library(Rcpp)
sourceCpp("approximation/count_distances_exdiag.cpp")
mcc <- 4


## query points x_i which are within distance r of x

sigma2 <- 0.05
sigma2_tr <- 0.05

p <- 10
n <- 1e5
mus <- randn(n, p)

muhs <- mus + sqrt(sigma2_tr) * randn(n, p)
ys <- mus + sqrt(sigma2) * randn(n, p)
rSqs <- rowSums((ys - muhs)^2)

t1 <- proc.time()
counts <- countDistEx(muhs, ys, rSqs)
(cpp_time <- proc.time() - t1)

1 - mean(counts != 0) ## accuracy

## takes 290s for 1e5

## naive way
# t1 <- proc.time()
# dd <- pdist2(muhs, ys)^2
# counts0 <- sapply(1:n, function(i) sum(dd[-i, i] < rSqs[i]))
# (naive_time <- proc.time() -t1)
# 
# sum(counts != counts0)
# 1 - mean(counts != 0) ## accuracy
# rbind(cpp_time, naive_time)
