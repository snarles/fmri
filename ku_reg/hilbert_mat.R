library(numbers)

hmat <- function(n) {
  ans <- matrix(0, n, n)
  1/(row(ans) + col(ans) - 1)
}

ihmat <- function(n) {
  ans <- matrix(0, n, n)
  #solve(1/(row(ans) + col(ans) - 1))
  i <- as.numeric(row(ans))
  j <- as.numeric(col(ans))
  matrix(
    (-1)^(i+j) * (i + j - 1) *
    choose(n + i - 1, n - j) *
      choose(n + j - 1, n - i) *
      choose(i + j - 2, i - 1)^2,
    n, n)
}


hmat(3)


sapply(1:10, function(i) sum(solve(hmat(i))))

sapply(1:20, function(i) sum(ihmat(i)))

## positive parts
sapply(1:4, function(i) sum(pmax(ihmat(i), 0)))
## negative parts
sapply(1:4, function(i) sum(pmin(ihmat(i), 0)))


sum(pmax(ihmat(4), 0))
sum(pmin(ihmat(4), 0))


sum(pmax(ihmat(5), 0))
primeFactors(478025)

sum(pmin(ihmat(5), 0))
