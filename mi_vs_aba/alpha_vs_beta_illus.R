## Get the convex set of all alpha, beta pairs


get_q_ints <- function(qs, k, naive = FALSE, mc.reps = 1e4) {
  if (naive) {
    pmf <- qs/length(qs)
    lala <- sapply(1:mc.reps, function(i) {
      xs <- sample(length(pmf), k, replace = TRUE)
      ca <- max(qs[xs])/k
      ci <- mean(qs[xs] * log(qs[xs]))
      avg <- mean(qs[xs])
      c(avg, ci, ca)
    })
    return(rowMeans(lala))
  }
  qs <- sort(qs)
  d <- length(qs)
  c(mean(qs), mean(qs * log(qs)),
    sum((1/k) * qs *
          (
            ((1:d)/d)^k - 
              ((0:(d-1))/d)^k
          )
        )
    )
}


## dimension of pmf
d <- 5; k <- 2
pmf <- runif(d); qs <- sort(pmf/sum(pmf) * d)
qs
get_q_ints(qs, k, TRUE)
get_q_ints(qs, k)
