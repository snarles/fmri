source("approximation/gaussian_identity.R")
source("approximation/large_deviations_source.R")

#### 
##  Does the mc_ident3 (large-sigma approximation) work for large K?
####

## sigma2 = 10 ##
p <- 10; sigma2 <- 10; L <- 5; mc.reps <- 1e5
mc_ident(p, sigma2, L, mc.reps)
mc_ident2(p, sigma2, L, mc.reps)
mc_ident3(p, sigma2, L, mc.reps)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE, use.qchisq=TRUE)

p <- 20; sigma2 <- 10; L <- 30; mc.reps <- 1e5
mc_ident3(p, sigma2, L, mc.reps)
mc_ident(p, sigma2, L, mc.reps)
mc_ident2(p, sigma2, L, mc.reps)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE, use.qchisq=TRUE)

p <- 80; sigma2 <- 10; L <- 200; mc.reps <- 1e5
mc_ident3(p, sigma2, L, mc.reps)
mc_ident2(p, sigma2, L, mc.reps)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE, use.qchisq=TRUE)

p <- 100; sigma2 <- 10; L <- 300; mc.reps <- 1e4
mc_ident3(p, sigma2, L, mc.reps)
mc_ident2(p, sigma2, L, mc.reps)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE)

## sigma2 = 100 ##
p <- 90; sigma2 <- 100; L <- 5; mc.reps <- 1e4
mc_ident(p, sigma2, L, mc.reps)
mc_ident2(p, sigma2, L, mc.reps)
mc_ident3(p, sigma2, L, mc.reps)
mc_ident4(p, sigma2, L, mc.reps, exact = TRUE)


####
##  Does L = 1/sigma2 ^(d/2) preserve rate??
####
p <- 5; L0 <- 10; sigma20 <- 0.1;
1 %>% {mc_ident3(p, sigma20/., floor(L0 * .^(d/2)) )}
