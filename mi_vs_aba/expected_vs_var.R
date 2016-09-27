### Q quantile point of view


## Check max formula 
# k <- 3
# qs <- sort(rbeta(100, 2, 3))
# sum(qs * seq(0, 1, length.out = length(qs))^(k-1))/length(qs)
# mean(sapply(1:1000, function(i) max(sample(qs, k, TRUE))))/k


get_q_ints <- function(qs, k, naive = FALSE) {
  qs <- sort(qs)
  if (naive) {
    qsamp <- matrix(sample(qs, 1e4 * k, replace = TRUE), 1e4, k)
    eQ <- mean(qsamp)
    ms <- apply(qsamp, 1, max)/k
    maxQ <- mean(ms)
    varQ <- var(ms)
    return(c(eQ, maxQ, varQ))
  }
  eQ <- sum(qs)/length(qs)
  maxQ <- sum(qs/k * k * seq(0, 1, length.out = length(qs))^(k-1))/length(qs)
  maxQ2 <- sum((qs/k)^2 * k * seq(0, 1, length.out = length(qs))^(k-1))/length(qs)
  varQ <- maxQ2 - (maxQ)^2
  c(eQ, maxQ, varQ)
}

normalize_qs <- function(qs) {
  sort(qs)/sum(qs) * length(qs)
}

qs_par <- function(c1, c2, k, l.out = 100) {
  xs <- seq(0, 1, length.out = l.out)
  qs <- pmax(c1 - c2/(xs ^ (k-1)), 0)
  qs
}

qs_par_c1 <- function(c1, k, l.out = 100) {
  c2 <- find_c2(c1, k)
  qs_par(c1, c2, k, l.out)
}

compute_eQ <- function(c1, c2, k) {
  if ((c2/c1) > 1) {
    return(0)
  }
  c1 * (1 - (c2/c1)^(1/(k-1))) - c2/(2 - k) * (1 - (c2/c1)^((2-k)/(k-1)))
}

find_c2 <- function(c1, k, lower = 1e-10) {
  ff <- function(c2) {
    compute_eQ(c1, c2, k) - 1
  }
  uniroot(ff, interval = c(lower, c1))$root
}

k <- 5
reso <- 1000
qs <- rbeta(reso, runif(1) * 10, runif(1) * 10)

pts <- do.call(rbind, lapply(seq(1.1, 100, by = 0.1), function(c1) {
  qs <- qs_par_c1(c1, k, l.out = 1e4)
#  plot(qs, type = "l")
  get_q_ints(qs, k)
}))
# 
# pts <- sapply(1:10000, function(i) {


ptsR <- do.call(rbind, lapply(seq(1.1, 100, by = 0.1), function(i) {
  qs <- rbeta(reso, runif(1) * 10, runif(1) * 10)
  qs <- normalize_qs(qs)
  #  plot(qs, type = "l")
  get_q_ints(qs, k)
}))
plot(ptsR[, 2], sqrt(ptsR[, 3]), pch = ".")
lines(pts[, 2], sqrt(pts[, 3]), type = "l", col = "red")

  
  #   qs <- normalize_qs(qs)
#   get_q_ints(qs, k)
# })
# 
# pts2 <- sapply(1000 * (1:10000/10000), function(x) {
#   qs <- qs_par(x, k, reso)
#   get_q_ints(qs, k)
# })
# 
# plot(t(pts)[, -1], pch= ".")
# lines(t(pts2)[, -1], col = "red")
# 
# plot(qs_par(d = 2, k = 3), type = "l")
# plot(qs_par(d = 2, k = 4), type = "l")
# plot(qs_par(d = 2, k = 4), type = "l")
# plot(qs_par(d = 10, k = 4), type = "l")

