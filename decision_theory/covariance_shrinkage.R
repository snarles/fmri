####
##  Covariance shrinkage under gaussian loss
####

gaussian_cov_loss <- function(Sigma, SigmaH)
  sum(diag(solve(SigmaH, Sigma))) + log(det(SigmaH))

####
##  For given sigma, best off-diagonal shrinkage param
####

n <- 20; p <- 30
Sigma <- cov(randn(40, p))
dim(Sigma)
mc.reps <- 10
lambdas <- seq(0.05, 1, by = .05)
losses <- matrix(0, mc.reps, length(lambdas))

for (i in 1:mc.reps) {
  X <- mvrnorm(n, rep(0, p), Sigma)
  SigmaH <- (t(X) %*% X)/n
  for (j in 1:length(lambdas)) {
    lambda <- lambdas[j]
    SigmaR <- lambda * diag(diag(SigmaH)) + (1-lambda) * SigmaH
    losses[i, j] <- gaussian_cov_loss(Sigma, SigmaR)
  }
}
plot(lambdas, colMeans(losses), type = "l")
