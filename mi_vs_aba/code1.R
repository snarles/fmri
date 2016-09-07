#####
##  Relationship between ABA and MI (discrete dist)
#####

colnorm <- function(mat, nits = 20) {
  for (i in 1:nits) {
    mat <- mat/rowSums(mat)/nrow(mat)
    mat <- t(t(mat)/colSums(mat))/ncol(mat)
  }
  mat
}

mi_mat <- function(mat) {
  px <- rowSums(mat)
  py <- colSums(mat)
  pxy <- t(t(px)) %*% t(py)
  as.numeric(sum(mat * log(mat/pxy)))
}

aba_mat_naive <- function(mat, k=2, nits = 1e3) {
  px <- rowSums(mat)
  trials <- sapply(1:nits, function(i) {
    inds <- sample(1:nrow(mat), k, replace = TRUE, prob = px)
    z <- sample(1:k, 1)
    y <- sample(1:ncol(mat), 1, prob = mat[inds[z], ]/sum(mat[inds[z], ]))
    (which.max(mat[inds, y]) == z)[1]
  })
  mean(trials)
}

library(AlgDesign)

aba_mat <- function(mat, k, ntr = 100) {
  p <- nrow(mat)
  # dat <- gen.factorial(levels = p, nVars = k, factors = "all")
  # desT <- optFederov(~., dat, nTrials = ntr)
  # dd <- desT$design
  # table(as.character(as.matrix(dd)))
  px <- rowSums(mat)
  trials <- sapply(1:nits, function(i) {
    inds <- sample(1:nrow(mat), k, replace = TRUE, prob = px)
    mat2 <- mat[inds, ]
    mat2 <- 1/k * mat2/rowSums(mat2)
    sum(apply(mat2, 2, max))
  })
  mean(trials)
}

library(pracma)
p <- 20
mat <- exp(randn(p)) + 10 * eye(p)
mat <- colnorm(mat, nits = 20)
rowSums(mat)

mi_mat(mat)
aba_mat(mat, k = 3, ntr = 1000)
aba_mat_naive(mat, k = 3, nits = 1e4)
