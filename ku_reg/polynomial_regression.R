## variance of polynomial regression

library(GenSA)

lagwts <- function(xs, x) {
  mat1 <- xs %*% t(rep(1, length(xs)))
  diff <- t(mat1) - mat1
  diff2 <- diff / (x - xs)
  diag(diff2) <- 1
  invw <- apply(diff2, 2, prod)
  1/invw
}

varlag <- function(xs, x) {
  sum(lagwts(xs, x)^2)
}

bestvar <- function(d, ulim) {
  ff <- function(x) { 
    varlag(x, 1)
  }
  res <- 
    GenSA(par = NULL, lower = rep(0, d+1), upper = rep(ulim, d+1), fn = ff)
  print(res$value)
  sort(res$par)
}

ppdf <- 5 ## number of points per df
n <- ppdf * (d + 1) ## number of points
d <- 5  ## degree
xs <- (1:n)/(n+1)
xmat <- sapply(0:d, function(k) xs^k)
xvec <- rep(1, d + 1)
wvec0 <- (t(xvec) %*% solve(t(xmat) %*% xmat, t(xmat)))[1, ]
(ols_var <- sum(wvec0^2))
(t(xvec) %*% solve(t(xmat) %*% xmat, xvec))[1]


# test of lagwts
# inds <- sample(1:n, d + 2)
# w <- lagwts(xs[inds[-1]], xs[inds[1]])
# sum(w * xs[inds[-1]]^d)
# xs[inds[1]]^d


####
## approximating variance with choose 1 from each partition
## this is bad!
####

mat <- AlgDesign::gen.factorial(rep(ppdf, d + 1), center = FALSE)
lala <- t(mat) + (0:d) * ppdf
ws_mat <- apply(lala, 2, function(v) {
  ans <- rep(0, n)
  ans[v] <- lagwts(xs[v], 1)
  ans
})
avw <- rowMeans(ws_mat)
sum(avw^2)
vars <- colSums(ws_mat^2)
sort(vars)[1:100]

####
##  approximating with disjoint subsets!
####

mat <- rep(ppdf, ppdf) %*% t(0:d)
lala <- t(mat + 1:ppdf)
lala
ws_mat <- apply(lala, 2, function(v) {
  ans <- rep(0, n)
  ans[v] <- lagwts(xs[v], 1)
  ans
})
avw <- rowMeans(ws_mat)
sum(avw^2)
vars <- colSums(ws_mat^2)
vars

####
## approximating variance with n choose (d + 1) sums
####

lala <- combn(1:n, d + 1)
ws_mat <- apply(lala, 2, function(v) {
  ans <- rep(0, n)
  ans[v] <- lagwts(xs[v], 1)
  ans
})
avw <- rowMeans(ws_mat)
sum(avw^2)
sum(ws_mat[, 2])

vars <- colSums(ws_mat^2)
hist(vars)
lala[, order(vars)[1:10]]
lala[, order(-vars)[1:10]]
sort(vars)[1:100]

ncombs <- 1:1000 * 10
varis <- sapply(ncombs, function(k) {
  avw2 <- rowMeans(ws_mat[, order(vars)[1:k], drop = FALSE])
  sum(avw2^2)
})

(d+1)^ppdf

plot(ncombs, varis, type = "l", ylim = c(0, 3 * ols_var))
abline(h = ols_var, col = "red")
