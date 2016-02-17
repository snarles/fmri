####
##  Bootstrap testing mu^2 = mu^2
##  One dimensional case!!
####

data_model <- function(muA, muB, sigma2A, sigma2B) {
  sampler <- function(nA, nB) {
    noise_mult <- 1
    if (nA + nB == 0) {
      noise_mult <- 0
      nA <- 2
      nB <- 2
    }
    xA <- rnorm(nA, mean = muA, sd = noise_mult * sqrt(sigma2A))
    xB <- rnorm(nB, mean = muB, sd = noise_mult * sqrt(sigma2B))
    dat <- rbind(cbind(0, xA), cbind(1, xB))
    dat
  }
  # force eval
  lala <- sampler(10, 10)
  sampler
}

stat.diff <- function(dat) {
  xA <- dat[dat[, 1]==0, 2]
  xB <- dat[dat[, 1]==1, 2]
  muA <- mean(xA); sigma2A <- var(xA)
  muB <- mean(xB); sigma2B <- var(xB)
  mu2A <- muA^2 - sigma2A/length(xA)
  mu2B <- muB^2 - sigma2B/length(xB)
  mu2A - mu2B
}

boot_sampler <- function(dat) {
  nX <- sum(dat[, 1] == 0); nY <- sum(dat[, 1] == 1);
  nX0 <- nX; nY0 <- nY
  rawX <- dat[dat[, 1]==0, -1]
  rawY <- dat[dat[, 1]==1, -1]
  sampler <- function(nX = nX0, nY = nY0) {
    if (nX == 0 & nY == 0) {
      return(dat)
    }
    newX <- rawX[sample(nX, nX, TRUE)]
    newY <- rawY[sample(nY, nY, TRUE)]
    newdat <- rbind(cbind(0, newX), cbind(1, newY))
    newdat
  }
  # force eval
  sampler()
  sampler
}

sampling_dist <- function(sampler, theta, nX, nY, mc.reps = 1e3, samples = FALSE) {
  # true params
  dat0 <- sampler(0, 0)
  theta0 <- theta(dat0)
  # get population
  thetas <- lapply(1:mc.reps, function(i) theta(sampler(nX, nY)))    
  thetas <- do.call(cbind, thetas)
  if (samples) {
    return(thetas)
  }
  diffs <- thetas - theta0
  bias <- rowMeans(diffs)
  sdv <- apply(diffs, 1, sd)
  skew <- rowSums(diffs^3)/(rowSums(diffs^2))^(3/2)
  list(bias = bias, sdv = sdv, skew = skew)
}

jackknife <- function(res, theta) {
  theta0 <- theta(res)
  n <- dim(res)[1]
  thetas <- list()
  for (i in 1:n) {
    dat2 <- res[-i, , drop = FALSE]
    thetas[[i]] <- theta(dat2)
  }
  thetas <- do.call(cbind, thetas)
  theta_bar <- rowMeans(thetas)
  diffs <- thetas - theta_bar
  bias <- (n-1) * (theta_bar - theta0)
  sdv <- sqrt((n-1) * rowMeans((thetas - theta_bar)^2))
  skew <- rowSums(diffs^3)/(rowSums(diffs^2))^(3/2)  
  list(bias = bias, sdv = sdv, skew = skew)
}

inv_alpha_bca <- function(alpha_bca, z0, acc) {
  if (alpha_bca == 1) {
    return(1)
  }
  za <- qnorm(alpha_bca)
  pnorm(
    ((1 - acc * z0) * (za - z0) - z0)/
      (1 + acc * (za - z0))
  )
}

inverse_bca_test <- function(res, theta, mc.reps = 1000) {
  jj <- jackknife(res, theta)
  acc <- jj$skew/6
  boot_s1 <- boot_sampler(res)
  nX <- sum(res[, 1] == 0); nY <- sum(res[, 1] == 1)
  thetas <- sampling_dist(boot_s1, theta, nX, nY, mc.reps, TRUE)
  theta0 <- theta(res)
  z0s <- sapply(1:length(theta0), function(i) {
    qnorm(sum(thetas[i, ] < theta0[i])/mc.reps)
  })
  pvs <- sapply(1:length(theta0), function(i) {
    v <- thetas[i, ]
    p1 <- (sum(v > 0) + 1)/(mc.reps + 1)
    p2 <- (sum(v < 0) + 1)/(mc.reps + 1)
    p1_mod <- inv_alpha_bca(p1, z0s[i], acc[i])
    p2_mod <- inv_alpha_bca(p2, -z0s[i], acc[i])
    2 * pmin(p1_mod, p2_mod)
  })
  min(pvs) * length(pvs)
}



####
##  Test out the functions
####

h0 <- data_model(0, 0, 1, 1)
h0v <- data_model(0,0,1, 10)
h1 <- data_model(0, 1, 1, 1)

stat.diff(h0(10, 10))
stat.diff(h0(1000, 1000))
stat.diff(h0(0,0))

stat.diff(h0v(10, 10))
stat.diff(h0v(1000, 1000))
stat.diff(h0v(0,0))

stat.diff(h1(10, 10))
stat.diff(h1(1000, 1000))
stat.diff(h1(0,0))

(dat <- h0(5, 5))
bs <- boot_sampler(dat)
bs(10, 10)

####
##  Investigate effect of bootstrap resampling
####

mc.reps <- 1000

h0 <- data_model(0, 0, 1, 1)
h0v <- data_model(0,0,1, 10)
h1 <- data_model(0, 1, 1, 1)

nX <- 100; nY <- nX

vals1 <- sampling_dist(h0, stat.diff, nX, nY, mc.reps, TRUE)
vals2 <- sampling_dist(boot_sampler(h0(nX, nY)), stat.diff, nX, nY, mc.reps, TRUE)
layout(t(t(1:2))); hist(vals1); hist(vals2)

vals1 <- sampling_dist(h0v, stat.diff, nX, nY, mc.reps, TRUE)
vals2 <- sampling_dist(boot_sampler(h0v(nX, nY)), stat.diff, nX, nY, mc.reps, TRUE)
layout(t(t(1:2))); hist(vals1); hist(vals2)

vals1 <- sampling_dist(h1, stat.diff, nX, nY, mc.reps, TRUE)
vals2 <- sampling_dist(boot_sampler(h1(nX, nY)), stat.diff, nX, nY, mc.reps, TRUE)
layout(t(t(1:2))); hist(vals1); hist(vals2)

####
##  How does inverse BCA do?
####
mc.reps = 1000
nX <- 100; nY <- nX 
res0 <- h0(nX, nY)
res0v <- h0v(nX, nY)
res1 <- h1(nX, nY)

c(inverse_bca_test(res0, stat.diff, mc.reps), 
  inverse_bca_test(res0v, stat.diff, mc.reps),
  inverse_bca_test(res1, stat.diff, mc.reps))

