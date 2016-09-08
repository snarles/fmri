### Q quantile point of view


## Check max formula 
# k <- 3
# qs <- sort(rbeta(100, 2, 3))
# sum(qs * seq(0, 1, length.out = length(qs))^(k-1))/length(qs)
# mean(sapply(1:1000, function(i) max(sample(qs, k, TRUE))))/k

get_q_ints <- function(qs) {
  qs <- sort(qs)
  c(sum(qs), sum(qs * log(qs)),
    sum(qs * seq(0, 1, length.out = length(qs))^(k-1)))/length(qs)
}

normalize_qs <- function(qs) {
  sort(qs)/sum(qs) * length(qs)
}

qs <- normalize_qs(rbeta(100, 1, 5))
get_q_ints(qs)

qs_par <- function(d, l.out = 100) {
  xs <- seq(0, 1, length.out = l.out)
  normalize_qs(exp( - d * xs))
}

pts <- sapply(1:1000, function(i) {
  qs <- rbeta(100, runif(1) * 5, runif(1) * 5)
  qs <- normalize_qs(qs)
  get_q_ints(qs)
})

pts2 <- sapply(10 * (1:1000/1000), function(x) {
  qs <- qs_par(x)
  get_q_ints(qs)
})

plot(t(pts)[, -1], pch= ".")
lines(t(pts2)[, -1], col = "red")

