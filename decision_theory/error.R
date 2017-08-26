####
##  Risk formula for Sigma = I incorrect
##  Charles Zheng
####


mI <- function(z, gamma) {
  -((z + gamma - 1) + sqrt((z + gamma-1)^2 - 4 * gamma * z))/
    (2 * gamma * z)
}
mId <- function(z, gamma) {
  -mI(z, gamma)/z + 1/(2 * gamma * z) *
    (-1 - (z- gamma - 1)/sqrt((z + gamma - 1)^2 - 4 * gamma*z))
}
risk <- function(alpha2, gamma, lambda) {
  1 + gamma * mI(-lambda, gamma) - lambda * (gamma - lambda * alpha2) * mId(-lambda, gamma)
}

## check that mId is correct
z <- -rexp(1); gamma <- rexp(1)
mId(z, gamma)
numDeriv::grad(function(z) mI(z, gamma), z)

## here is the problem: risk is not minimized by lambda = gamma/alpha2

gamma <- rexp(1); alpha2 <- rexp(1)
risk(alpha2, gamma, gamma/alpha2)
risk(alpha2, gamma, 2*gamma/alpha2)

## a plot
lambdas <- gamma/alpha2 * seq(.1, 10, by = .05)
plot(sqrt(lambdas), sapply(lambdas, function(x) risk(alpha2, gamma, x)), type = "l")
abline(v = sqrt(gamma/alpha2))
