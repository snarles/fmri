library(AlgDesign)
source("info_theory_sims/sim3source.R")

str2vec <- function(s) {
  as.numeric(substring(s, seq(1,nchar(s),1), seq(1,nchar(s),1)))
}

eta <- function(x) log(1 + exp(x))

pr_grid <- function(y, Bmat, reso = 5) {
  Btilde <- Bmat %*% diag(-y)
  xs <- as.matrix(gen.factorial(rep(reso, p))/reso * 2 * sqrt(2 * log(reso)))
  delta <- xs[2, 1] - xs[1, 1]
  mus <- xs %*% Btilde
  nms <- rowSums(xs^2)/2+ rowSums(eta(mus))
  delta^p * (1/sqrt(2*pi))^p * dim(xs)[1] * meanexp(-nms)
}

p <- 2
q <- 3
Bmat <- randn(p, q)
mc.reps <- 1e6
X <- randn(mc.reps, p)
ps <- 1/(1 + exp(-X %*% Bmat))
Y <- (rand(mc.reps, q) < ps) + 0
ylabs <- apply(Y, 1, function(v) paste(v, collapse = ""))
tab <- table(ylabs)
tab <- tab/sum(tab)

#(s <- "111")
(s <- sample(ylabs, 1))
y <- 2 * str2vec(s) - 1
(pr_emp <- tab[s])
(pr_true <- pr_grid(y, Bmat, 300))
