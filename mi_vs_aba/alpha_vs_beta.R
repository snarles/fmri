### Q quantile point of view


## Check max formula 
# k <- 3
# qs <- sort(rbeta(100, 2, 3))
# sum(qs * seq(0, 1, length.out = length(qs))^(k-1))/length(qs)
# mean(sapply(1:1000, function(i) max(sample(qs, k, TRUE))))/k


get_q_ints <- function(qs, k) {
  qs <- sort(qs)
  c(sum(qs), sum(qs * log(qs)),
    sum(qs * seq(0, 1, length.out = length(qs))^(k-1)))/length(qs)
}

normalize_qs <- function(qs) {
  sort(qs)/sum(qs) * length(qs)
}

qs_par <- function(d, k, l.out = 100) {
  xs <- seq(0, 1, length.out = l.out)
  normalize_qs(exp( d * xs^(k-1)))
}

aba_par <- function(d, k, l.out, inds = 3) {
  xs <- seq(0, 1, length.out = l.out)
  get_q_ints(normalize_qs(exp( d * xs^(k-1))), k)[inds]
}

find_par_aba <- function(aba, k, l.out = 1000, init = 2, nits = 20, info = TRUE) {
  lb <- 0
  ub <- init
  while(aba_par(ub, k, l.out) < aba) {
    ub <- ub * 2
  }
  for (i in 1:nits) {
    d <- (lb + ub)/2
    if (aba_par(d, k, l.out) < aba) {
      lb <- d
    }
    if (aba_par(d, k, l.out) > aba) {
      ub <- d
    }
  }
  if (info) {
    d <- (lb + ub)/2
    return(list(d = d, res = aba_par(d, k, l.out, 2:3)))
  }
  (lb + ub)/2
}


# k <- 5
# reso <- 1000
# 
# pts <- sapply(1:10000, function(i) {
#   qs <- rbeta(reso, runif(1) * 10, runif(1) * 10)
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

