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

induction_mult <- function(n) {
  ans <- matrix(0, n, n)
  #solve(1/(row(ans) + col(ans) - 1))
  i <- as.numeric(row(ans))
  j <- as.numeric(col(ans))
  matrix((n+j)/(n-i+1) * (n+i)/(n-j+1) - 1,n,n)
}

ih_induction_mult <- function(n) {
  ans <- matrix(0, n, n)
  #solve(1/(row(ans) + col(ans) - 1))
  i <- as.numeric(row(ans))
  j <- as.numeric(col(ans))
  matrix((-1)^(i+j) *
           factorial(i+n-1) *
           factorial(j+n-1)/
           factorial(i-1)^2/
           factorial(j-1)^2/
           factorial(n-i+1)/
           factorial(n-j+1),n,n)
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

# final rows
sapply(1:10, function(n) sum(ihmat(n)[n, ]))
sapply(1:10, function(n) (2*n-1)* choose(2*n-2, n-1))

# nn element
sapply(1:10, function(n) ihmat(n)[n, n])
sapply(1:10, function(n) (2*n-1)* choose(2*n-2, n-1)^2)

## more identities
sapply(1:10, function(n) sum(ihmat(n) * induction_mult(n)))
sapply(1:10, function(n) sum(ihmat(n) * induction_mult(n)))

list(ihmat(5) * induction_mult(5),
     ihmat(6))
list(ihmat(5) * induction_mult(5),
     (2*5 + 1) * ih_induction_mult(5))
sapply(1:5, function(n) sum(ih_induction_mult(n)))
sapply(1:5, function(n) 1+choose(2*n, n)*(choose(2*n,n) - 2))

## looking for patterns
sapply(1:5, function(n) sum(pmax(ih_induction_mult(n), 0)))
sapply(1:5, function(n) sum(pmin(ih_induction_mult(n), 0)))

