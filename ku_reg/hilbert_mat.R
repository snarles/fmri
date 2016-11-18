hmat <- function(n) {
  ans <- matrix(0, n, n)
  1/(row(ans) + col(ans) - 1)
}

hmat(3)


nseq <- sapply(1:10, function(i) sum(solve(hmat(i))))
nseq
