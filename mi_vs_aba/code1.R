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

## E[Z(Y) log (Z(Y))]
mi_mat2 <- function(mat) {
  mat <- mat * ncol(mat) * nrow(mat)
  mean(mat * log(mat))
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

aba_mat_trials <- function(mat, k, ntr = 100, rtrials = TRUE) {
  p <- nrow(mat)
  # dat <- gen.factorial(levels = p, nVars = k, factors = "all")
  # desT <- optFederov(~., dat, nTrials = ntr)
  # dd <- desT$design
  # table(as.character(as.matrix(dd)))
  px <- rowSums(mat)
  trials <- sapply(1:ntr, function(i) {
    inds <- sample(1:nrow(mat), k, replace = TRUE, prob = px)
    mat2 <- mat[inds, ]
    mat2 <- 1/k * mat2/rowSums(mat2)
    sum(apply(mat2, 2, max))
  })
  if (rtrials) {
    return(trials)
  }
  mean(trials)
}

## (1/k) E[max_i Z_i(Y)]
aba_mat <- function(mat, k) {
  ## convert to continuous dist
  mat2 <- mat * nrow(mat) * ncol(mat)
  emaxs <- sapply(1:ncol(mat2), function(i) {
    vals <- mat2[, i]
    distr <- data.frame(support = sort(vals), pmf = 1/length(vals))
    expected_max(distr, k)
  })
  mean(emaxs)/k
}

pmf_max_k <- function(distr, k) {
  data.frame(support = distr$support, pmf = diff(c(0, cumsum(distr$pmf)^k)))
}

## takes a distribution (a data frame with support and pmf) and compute E[max_k]
expected_max <- function(distr, k, naive = FALSE) {
  if (naive) {
    lala <- sapply(1:10000, function(i) {
      bleh <- sample(distr$support, k, TRUE, distr$pmf)
      max(bleh)
    })
    return(mean(lala))
  }
  distr_k <- pmf_max_k(distr, k)
  tmax <- max(distr_k$support)
  tmax - sum(distr_k$pmf * (tmax - distr_k$support))
}



