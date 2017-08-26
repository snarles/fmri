####
##  General non-central chi squared
####

library(pracma); library(magrittr)

f2 <- function(x, y=0) sum((x-y)^2)

## probability that x in S^(p-1) has x[1] > 1-2*u
lcap_prob <- function(p, u) {
  log(1/2) + pbeta((1 - 2*u)^2, 1/2, (p+1)/2, lower.tail=FALSE, log.p = TRUE)
}

rgchisq0 <- function(n, Sigma, mu) {
  p <- dim(Sigma)[1]
  L <- chol(Sigma)
  x <- L %*% (randn(p, n) + mu)
  colSums(x^2)
}

rgchisq <- function(n, Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  x <- sqrt(ls) * (randn(p, n) + nu)
  colSums(x^2)
}

mgf_gchisq <- function(tt, Sigma, mu, log = FALSE) {
  p <- dim(Sigma)[1]
  tt <- t(t(tt))
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  tM <- repmat(tt, 1, p)
  lM <- repmat(t(ls), length(tt), 1)
  nM <- repmat(t(nu), length(tt), 1)
  if (log) {
    temp <- -1/2 * log(1 - 2*tM*lM) + (nM^2*lM*tM/(1 - 2*tM*lM))
    return(rowSums(temp))
  } else {
    temp <- 1/sqrt(1 - 2*tM*lM)*exp(nM^2*lM*tM/(1 - 2*tM*lM))
    return(apply(temp, 1, prod))
  }
}

## gives function handles for psi = log mgf, (d/dt) psi, (d2/dt^2) psi, (d psi)^{-1}
log_mgf_gchisq_ <- function(Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  psi <- function(tt)
    sum(-1/2 * log(1 - 2*tt*ls) + (nu^2*ls*tt/(1 - 2*tt*ls)))
  dpsi <- function(tt) {
    ##numDeriv::grad(psi, tt)
    sum(ls*(1 + nu^2)/(1-2*ls*tt)) + sum(2*nu^2*ls^2*tt/(1 - 2*ls*tt)^2)
  }
  d2psi <- function(tt) {
    ##numDeriv::hessian(psi, tt)[1]
    sum(2*ls^2*(1+2*nu^2)/(1-2*ls*tt)^2 + 8*ls^3*nu^2*tt/(1-2*ls*tt)^3)
  }
  dpsi_inv <- function(x, ...) {
    ff <- function(tt) (dpsi(tt) - x)^2
    res <- optimize(ff, ...)
    res$minimum
  }
  list(psi = psi, dpsi = dpsi, d2psi = d2psi, dpsi_inv = dpsi_inv)
}


## fourier inversion
finv <- function(fs, cs, x) {
  diffs <- Im(fs)[-1] - Im(fs)[-length(fs)]
  fM <- repmat(t(fs), length(x), 1)
  cM <- repmat(t(cs), length(x), 1)
  xM <- repmat(t(t(x)), 1, length(fs))
  vals <- exp(-fM * xM) * cM
  vals2 <- (vals[, -1, drop = FALSE] 
            + vals[, -length(fs), drop = FALSE])/2
  Re(rowSums(vals2 * diffs))/(2 * pi)
}

####
##  Upper bound using Markov
####

## exponential bound on Pr[<x]
mb_gchisq <- function(tt, x, Sigma, mu) {
  ms <- exp(mgf_gchisq(tt, Sigma, mu, TRUE))
  ms/exp(tt * x)
}

## derivative of log. expo. bound wrt t
prox_dleb <- function(tt, x, Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  u <- 1/tt
  # actual
  sum(ls*(1 + nu^2)/(1-2*ls*tt)) + sum(2*nu^2*ls^2*tt/(1 - 2*ls*tt)^2) - x
  # second-order
  (-p/2)*u + sum((1+2*nu+3*nu^2)/ls)/4*u^2 - x
  # first-order
  (-p/2) * u - x
}

## exponential bound on log Pr[<x] plugging in t = -p/2x
# lmb_gchisq <- function(x, Sigma, mu) {
#   log(mb_gchisq(-.5*dim(Sigma)[1]/x, x, Sigma, mu))
# }

lmb_gchisq <- function(x, Sigma, mu) {
    p <- dim(Sigma)[1]
    res <- eigen(Sigma); ls <- res$values
    nu <- as.numeric(t(res$vectors) %*% mu)
    -1/2 * sum(log(1 + ls*p/x)) + p/2 * 
      (1 - sum(nu^2*ls/(1 + ls*p/x))/x)
}

qlmb_gchisq <- function(lprob, Sigma, mu, intv = c(1e-10, 1e3)) {
  ff <- function(x) (lmb_gchisq(x, Sigma, mu) - lprob)^2
  res <- optimise(ff, interval=intv)
  res$minimum
}

####
##  Exponential tilting
####

quadfmla <- function(a, b, c, s = c(-1, 1), pos = FALSE) {
  disc <- b^2 - 4*a*c
  if (pos) { disc <- pmax(0, disc) }
  (-b + s*sqrt(disc))/(2 * a)  
}
lse <- function(x) max(x) + log(sum(exp(x - max(x))))

reso <- 1000
x <- fms$dpsi(0) * .01

l_gchisq_tilt_ <- function(Sigma, mu, reso = 1e4, lb = -1e4) {
  p <- dim(Sigma)[1]
  fms <- log_mgf_gchisq_(Sigma, mu)
  ##exp(fms$psi(tt)) * dlta * sum(exp(-tt * xs) * sig * dchisq(xs * sig, df, ncp))
  lgchisq <- function(x) {
    tt <- fms$dpsi_inv(x, interval = c(lb, 0))
    mm <- fms$dpsi(tt); vv <- fms$d2psi(tt)
    (df <- min(p, 2*mm^2/(vv)))
    (sig <- quadfmla(-2*df, 4 * mm, -vv, 1, TRUE))
    (ncp <- pmax(0, mm/sig - df))
    xs <- x * (1:reso)/reso
    dlta <- xs[2] - xs[1]
    ds <- dchisq(xs / sig, df, ncp, TRUE)
    if (sum(is.na(ds)) > 0) {
      print(list(df = df, sig = sig, ncp = ncp, x = x))
      ds_e <<- ds
    }
    fms$psi(tt) + log(dlta) +
      lse(-tt * xs - log(sig) + ds)    
  }
  qgchisq <- function(lprob, lb=-4) {
    ff <- function(x) (lgchisq(exp(x)) - lprob)^2
    res <- optimise(ff, interval=c(lb, log(fms$dpsi(0) * 10)))
    exp(res$minimum)
  }
  list(lgchisq = lgchisq, qgchisq = qgchisq)
}

####
##  Lower bound using cap
####

log_vsph <- function(d, r=1)
  log(2/d) + (d/2)*log(pi) - lgamma(d/2) + d*log(r)
#c(exp(log_vsph(3)), 4/3 * pi)

## bound on log Pr[<x], creates function handle
cap_lb_ <- function(x, Sigma, mu) {
  p <- dim(Sigma)[1]
  res <- eigen(Sigma); ls <- res$values
  nu <- as.numeric(t(res$vectors) %*% mu)
  ## Compute height constants, etc
  nun <- nu/sqrt(f2(nu))
  nSn <- sum(nun^2/ls)
  #c(nSn, t(mu) %*% solve(Sigma, mu)/f2(mu))
  zh <- sqrt(x/nSn) * (nun/ls)
  h <- sqrt(x * nSn)
  ## Volume of ellipse
  vc <- log_vsph(p) + (p/2)*log(x) - (1/2)*sum(log(ls))
  #c(h, sum(zh * nun))
  d0 <- sqrt(f2(h * nun, zh))
  ## Function handle
  ff <- function(u) {
    ## volume of cap
    cv <- vc + lcap_prob(p, u)
    ml <- (sqrt(f2(nu)) + h*(1-2*u))^2 +
      (d0 + sqrt(x/min(ls))*sqrt(1-(1-2*u)^2))^2
    dens <- -(p/2) * log(2*pi) - (ml^2)/2
    dens + cv
  }
  ff
}

# ####
# ##  Check cap volume and density
# ####
# 
# p <- 3
# mu <- 2 * rnorm(p)
# Sigma <- 0*cov(randn(3*p, p)) + eye(p)
# s1 <- rgchisq(1e6, Sigma, mu)
# x <- s1[100]
# res <- eigen(Sigma); ls <- res$values
# nu <- as.numeric(t(res$vectors) %*% mu)
# cube_l_min <- 2 * sqrt(x/min(ls))
# 
# mc.s <- 1e5
# cube_l <- 1.1 * cube_l_min
# cube_vol <- cube_l^p
# pts <- cube_l * matrix(runif(mc.s * p) - .5, mc.s, p)
# 
# ## check volume
# mc.s <- 1e7
# pts <- cube_l * matrix(runif(mc.s * p) - .5, mc.s, p)
# #nms <- apply(pts, 1, function(v) t(v) %*% Sigma %*% v)
# nms <- rowSums((pts %*% sqrtm(Sigma))^2)
# #  # empirical volume
# mean(nms < x) * cube_vol
# #  # computed volume
# (vol_e <- exp(log_vsph(p) - .5 * sum(log(ls)) + (p/2) * log(x)))
# ee <- pts[nms < x, ]
# 
# ## Check cap height
# nun <- nu/sqrt(f2(nu))
# nSn <- sum(nun^2/ls)
# zh <- sqrt(x/nSn) * (nun/ls)
# h <- sqrt(x * nSn)
# mun <- mu/sqrt(f2(mu))
# hs <- ee %*% mun
# c(max(hs), h)
# d0 <- sqrt(f2(h * nun, zh))
# 
# ## Check cap vol
# u <- runif(1)/3
# #  #  empirical volume
# mean(hs > (1-2*u)*h) * vol_e
# #  #  computed volume
# exp(log_vsph(p) - .5 * sum(log(ls)) + (p/2) * log(x) + lcap_prob(p, u))
# 
# ## Check cap distance
# diffs <- t(-t(ee) + mu)
# par_dists <- (diffs %*% mun)^2
# dists <- rowSums(t(t(ee) - mu)^2)
# orth_dists <- dists - par_dists
# (maxdist <- max(dists[hs > (1-2*u)*h]))
# (sqrt(f2(nu)) - h*(1-2*u))^2 +  (d0 + sqrt(x/min(ls))*sqrt(1-(1-2*u)^2))^2
# c(max(par_dists[hs > (1-2*u)*h]), (sqrt(f2(nu)) - h*(1-2*u))^2)
# c(max(orth_dists[hs > (1-2*u)*h]), (d0 + sqrt(x/min(ls))*sqrt(1-(1-2*u)^2))^2)

# ####
# ##  Tests
# ####
# 
# p <- 10
# mu <- rnorm(p)
# Sigma <- cov(randn(2*p, p))
# #Sigma <- eye(p)
# #s1 <- rgchisq0(1e6, Sigma, mu)
# s1 <- sort(rgchisq(1e6, Sigma, mu))
# cdf <- function(x) sum(s1 < x)/length(s1)
# cap_par <- 0.01
# lf <- l_gchisq_tilt_(Sigma, mu)$lgchisq
# qf <- l_gchisq_tilt_(Sigma, mu)$qgchisq
# lf(qf(-200, -40))
# s1[50] %>% {c(cap_lb_(., Sigma, mu)(cap_par), log(cdf(.)), lmb_gchisq(., Sigma, mu))}
# s1[50] %>% {c(log(cdf(.)), lmb_gchisq(., Sigma, mu), lf(.))}
# 
# #pchisq(s1[50], p, f2(mu), log.p=TRUE) ## if Sigma==eye(p)
# 
# x <- s1[50]
# ff <- cap_lb_(x, Sigma, mu)
# (1:1e3/1e7) %>% plot(., sapply(., ff), type = "l")
