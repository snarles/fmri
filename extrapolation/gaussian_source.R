library(pracma)
library(class)

get_bayes_acc <- function(mus, mc.res = 10) {
  k <- nrow(mus)
  p <- ncol(mus)
  zs <- rep(1:k, each = mc.res)
  mus2 <- mus[zs, , drop = FALSE]
  ys <- mus2 + randn(k * mc.res, p)
  cl <- knn(train = mus, test = ys, cl = 1:k)
  sum(cl == zs)/length(zs)
}

get_vs <- function(mus, r) {
  k <- nrow(mus)
  p <- ncol(mus)
  zs <- rep(1:k, each = r)
  mus2 <- mus[zs, , drop = FALSE]
  ys <- mus2 + randn(k * r, p)
  dd <- -pracma::pdist2(mus, ys)
  rankconv <- apply(dd, 2, rank)
  Vs <- list()
  for (i in 1:k) {
    Vs[[i]] <- rankconv[i, 1:r + (r * (i-1))]  
  }
  Vs <- do.call(c, Vs)
  Vs
}

compute_g0 <- function(p, g.res, mc.res) {
  psM <- (1:mc.res - 0.5)/mc.res
  us <- (1:g.res - 0.5)/g.res
  ubreaks <- (0:g.res)/g.res
  nms2 <- qchisq(psM, df = p)
  ans <- 0 * us
  for (i in 1:mc.res) {
    nm <- nms2[i]
    qbreaks <- qchisq(1-ubreaks, p, ncp = nm)
    pbreaks <- pchisq(qbreaks/4, p)
    pcont <- pbreaks[-(g.res + 1)] - pbreaks[-1]
    ans <- ans + pcont/mc.res * g.res
  }
  ans
}





####
##  Testing code
####
p <- 1
k <- 10
mc.res <- 1000
g.res <- 100
r <- 2
mus <- randn(k, p)

t1 <- proc.time()
mus <- randn(5000, 1)
vs <- get_vs(mus, 1)
g.res <- 100
us <- (1:g.res - 0.5)/g.res
ubreaks <- (0:g.res)/g.res
res <- hist(vs/nrow(mus), breaks = ubreaks, plot = FALSE)
gu.emp <- res$counts/length(vs) * g.res
plot(us, gu.emp, type = "l")
proc.time() - t1

t1 <- proc.time()
gu0 <- compute_g0(1, 100, 1000)
proc.time() - t1
sum(gu0)
lines(us, gu0, col = "red")
