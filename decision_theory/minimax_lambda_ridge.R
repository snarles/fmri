source("decision_theory//ridge_error.R")
source("utils//avger.R")

mI <- function(z, gamma) {
  -((z + gamma - 1) + sqrt((z + gamma-1)^2 - 4 * gamma * z))/
    (2 * gamma * z)
}
mId <- function(z, gamma) {
  -mI(z, gamma)/z + 1/(2 * gamma * z) *
    (-1 - (z- gamma - 1)/sqrt((z + gamma - 1)^2 - 4 * gamma*z))
}


z <- -rexp(1); gamma <- rexp(1)
mId(z, gamma)
numDeriv::grad(function(z) mI(z, gamma), z)

risk <- function(alpha2, gamma, lambda) {
  1 + gamma * mI(-lambda, gamma) + lambda * (gamma - lambda * alpha2) * mId(-lambda, gamma)
}
risk_ <- function(alpha2, gamma) {
  ff <- function(lambda) risk(alpha2, gamma, lambda)
  ff
}
opt_risk <- function(alpha2, gamma, naive = FALSE) {
  if (naive) {
    ff <- function(lambda) risk(alpha2, gamma, lambda)
    res <- optimise(ff, c(1e-5,1e5))
    return(res$objective)
  }
  .5 * (1 + (gamma-1)/gamma * alpha2 + sqrt((1-(gamma-1)/gamma*alpha2)^2 + 4*alpha2))
}

alpha2 <- rexp(1); gamma <- rexp(1)
risk(alpha2, gamma, gamma/alpha2)
avger(ridge_error_random_I, 100, 7, TRUE, alpha2, gamma, gamma/alpha2, naive = TRUE)
ridge_error_random_I(alpha2, gamma, gamma/alpha2, naive = TRUE)
risk(alpha2, gamma, 2 * gamma/alpha2)


