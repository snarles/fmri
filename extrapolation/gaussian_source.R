library(pracma)
library(class)

get_bayes_acc <- function(mus, sigma = 1, mc.res = 10) {
  k <- nrow(mus)
  p <- ncol(mus)
  zs <- rep(1:k, each = mc.res)
  mus2 <- mus[zs, , drop = FALSE]
  ys <- mus2 + sigma * randn(k * mc.res, p)
  cl <- knn(train = mus, test = ys, cl = 1:k)
  sum(cl == zs)/length(zs)
}

get_vs <- function(mus, sigma = 1, r = 1) {
  k <- nrow(mus)
  p <- ncol(mus)
  zs <- rep(1:k, each = r)
  mus2 <- mus[zs, , drop = FALSE]
  ys <- mus2 + sigma * randn(k * r, p)
  dd <- -pracma::pdist2(mus, ys)
  rankconv <- apply(dd, 2, rank)
  Vs <- list()
  for (i in 1:k) {
    Vs[[i]] <- rankconv[i, 1:r + (r * (i-1))]  
  }
  Vs <- do.call(c, Vs)
  Vs
}

get_us <- function(mus, sigma = 1, r = 1) {
  k <- nrow(mus)
  p <- ncol(mus)
  zs <- rep(1:k, each = r)
  mus2 <- mus[zs, , drop = FALSE]
  errs <- sigma * randn(k * r, p)
  ys <- mus2 + errs
  nms <- rowSums(ys^2)
  enms <- rowSums(errs^2)
  1 - pchisq(enms, p, ncp = nms)
}

compute_gu <- function(p, sigma = 1, g.res = 100, mc.res = 1e6) {
  us <- (1:g.res - 0.5)/g.res
  ubreaks <- (0:g.res)/g.res
  mus <- randn(mc.res, p)
  usamp <- get_us(mus, sigma, 1)
  res <- hist(usamp, breaks = ubreaks, plot = FALSE)
  gu.emp2 <- res$counts/length(usamp) * g.res
  list(gu = gu.emp2, mids = us, breaks = ubreaks)
}


####
##  Testing code
####
# p <- 1
# k <- 10
# mc.res <- 1000
# g.res <- 100
# r <- 2
# mus <- randn(k, p)
# 
# p <- 2
# t1 <- proc.time()
# mus <- randn(2000, p)
# vs <- get_vs(mus, 1)
# g.res <- 100
# us <- (1:g.res - 0.5)/g.res
# ubreaks <- (0:g.res)/g.res
# res <- hist(vs/nrow(mus), breaks = ubreaks, plot = FALSE)
# gu.emp <- res$counts/length(vs) * g.res
# plot(us, gu.emp, type = "l")
# proc.time() - t1
# 
# t1 <- proc.time()
# mus <- randn(1e6, p)
# usamp <- get_us(mus, 1)
# res <- hist(usamp, breaks = ubreaks, plot = FALSE)
# gu.emp2 <- res$counts/length(usamp) * g.res
# lines(us, gu.emp2, col = "red")
# proc.time() - t1
# 
# 
# sum(gu0)
# lines(us, gu0, col = "red")
# 
# res <- compute_gu(1, 0.1)
# plot(res$mids, res$gu, type = "l")
