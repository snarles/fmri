source("approximation/gaussian_identity.R")
source("approximation/gaussian_lc_source.R")
source("approximation/large_deviations_source.R")

#### 
##  Check fixed-L results
####

noise_const <- 2
L <- 5
mcK_I(cc=1/noise_const, Ks=L, mc.reps=1e4)

p <- 10; mc_ident(p, noise_const * p, L, 1e4)
p <- 20; mc_ident(p, noise_const * p, L, 1e4)
p <- 30; mc_ident(p, noise_const * p, L, 1e4)

p <- 10; mc_ident2(p, noise_const * p, L, 1e4)
p <- 20; mc_ident2(p, noise_const * p, L, 1e4)
p <- 30; mc_ident2(p, noise_const * p, L, 1e4)

p <- 10; mc_ident4(p, noise_const * p, L, 1e4, exact = TRUE)
p <- 10; mc_ident5(p, noise_const * p, L, 1e2)
p <- 10; mc_ident5(p, noise_const * p, L, 1e2, pchisq = pchisq_laplace)


noise_const <- .2
L <- 1e3
mcK_I(cc=1/noise_const, Ks=L, mc.reps=1e4)

p <- 1e3; mc_ident2(p, noise_const * p, L, 1e4)
mc_ident4(p, noise_const * p, L, 1e2)
mc_ident5(p, noise_const * p, L, 1e2)
mc_ident5(p, noise_const * p, L, 1e2, pchisq = pchisq_laplace)

p <- 2e3; mc_ident2(p, noise_const * p, L, 1e4)
mc_ident4(p, noise_const * p, L, 1e2)
mc_ident5(p, noise_const * p, L, 1e2)

p <- 3e3; mc_ident2(p, noise_const * p, L, 1e4)
mc_ident4(p, noise_const * p, L, 1)
mc_ident5(p, noise_const * p, L, 1)
mc_ident5(p, noise_const * p, L, 1, pchisq = pchisq_laplace)





noise_const <- .05
L <- 1e4
mcK_I(cc=1/noise_const, Ks=L, mc.reps=1e4)

p <- 1e3; mc_ident2(p, noise_const * p, L, 1e3)
mc_ident4(p, noise_const * p, L, 1e3)
mc_ident5(p, noise_const * p, L, 1e3)
mc_ident5(p, noise_const * p, L, 1, pchisq = pchisq_laplace)

p <- 2e3; mc_ident2(p, noise_const * p, L, 1e4)
mc_ident4(p, noise_const * p, L, 1e2)
mc_ident5(p, noise_const * p, L, 1e2)

p <- 3e3; mc_ident2(p, noise_const * p, L, 1e4)
mc_ident4(p, noise_const * p, L, 10)
mc_ident5(p, noise_const * p, L, 1e2)


# sigma2 <- noise_const * p
# y2 <- (1 + sigma2) * rchisq(1, df = p)
# K <- L; df <- p; ncp <- y2; nits <- 20

noise_const <- .0525
L <- 1e5; mc.reps <- 1e2
mcK_I(cc=1/noise_const, Ks=L, mc.reps=1e4)

p <- 1e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 2e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 4e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)

noise_const <- .0425
L <- 1e6; mc.reps <- 1e2
mcK_I(cc=1/noise_const, Ks=L, mc.reps=1e4)

p <- 1e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 2e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 4e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 8e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)

noise_const <- .0355
L <- 1e7; mc.reps <- 1e2
mcK_I(cc=1/noise_const, Ks=L, mc.reps=1e4)

p <- 1e2; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 1e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 2e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)

noise_const <- .02
L <- 1e12; mc.reps <- 1e3
mcK_I(cc=1/noise_const, Ks=L, mc.reps=1e4)
p <- 1e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)
p <- 2e3; mc_ident5(p, noise_const * p, L, mc.reps)
mc_ident5(p, noise_const * p, L, mc.reps, pchisq=pchisq_laplace)


## Fixed-sigma results: dimensionality-shortage effect?

## Formula is accurate for large sigma2 relative to ... dimensionality?
mc.reps <- 1e3
sigma2 <- 100; lc <- 0.01
p <- 1e2; 
mc_ident2(p, sigma2, floor(exp(lc * p)), mc.reps)

mc_ident5(p, sigma2, floor(exp(lc * p)), mc.reps)
mcK_I(cc=p/sigma2, Ks=floor(exp(lc * p)), mc.reps=1e4)
p <- 2e2; mc_ident5(p, sigma2, floor(exp(lc * p)), mc.reps)
mcK_I(cc=p/sigma2, Ks=floor(exp(lc * p)), mc.reps=1e4)
p <- 3e2; mc_ident5(p, sigma2, floor(exp(lc * p)), mc.reps)
mcK_I(cc=p/sigma2, Ks=floor(exp(lc * p)), mc.reps=1e4)
p <- 1e3; mc_ident5(p, sigma2, floor(exp(lc * p)), mc.reps)
mcK_I(cc=p/sigma2, Ks=floor(exp(lc * p)), mc.reps=1e4)

sigma2 <- 100; p <- 1e3; K<- 1e3
mcK_I(cc=p/sigma2, Ks=K, mc.reps=1e4)
mc_ident5(p, sigma2, K, mc.reps)

#mc_ident5(p, sigma2, K, mc.reps, pchisq_f=pchisq_laplace)
sigma2 <- 100; p <- 2e3; K<- 1e3
mcK_I(cc=p/sigma2, Ks=K, mc.reps=1e4)
mc_ident5(p, sigma2, K, mc.reps)



# Dbugging
sigma2 <- noise_const * p
K <- L; df <- p; nits <- 20
mc.reps <- 5
quants <- (1:mc.reps - .5)/mc.reps
y2s <- (1 + sigma2) * qchisq(quants, df = p)
## dbugging inner
K <- L-1
ncp <- y2s[1]
ux <- df + ncp; lx <- 0
const <- -log(2)
for (i in 1:nits) {
  (xc <- (ux + lx)/2)
  pchisq(xc, df, ncp)
  (val <- K * log_1_plus(-pchisq(xc, df, ncp)))
  if (val > const) {
    lx <- xc
  } else {
    ux <- xc
  }
  #print(xc)
}

## Error~~!
df <- 8000; ncp <- 2672870
pchisq(df + ncp, df, ncp)
pchisq(2010653, 8000, 2672870)
pchisq_laplace(2010653, 8000, 2672870)

