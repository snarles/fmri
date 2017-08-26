## testing RNG methods in R

master.seed <- 9999
region.seed <- 123
set.seed(master.seed * 1e5 + region.seed)
rnorm(10)
region.seed <- 124
set.seed(master.seed * 1e5 + region.seed)
rnorm(10)
