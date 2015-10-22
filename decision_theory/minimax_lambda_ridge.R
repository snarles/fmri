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
  1 + gamma * mI(-lambda, gamma) - lambda * (gamma - lambda * alpha2) * mId(-lambda, gamma)
}
risk2 <- function(alpha2, gamma, lambda) {
  Q <- (gamma - lambda - 1)^2 + 4 * gamma * lambda
  1/2 + alpha2/2 - alpha2/(2 * gamma) + alpha2*sqrt(Q)/(2*gamma) + 
    (gamma-lambda*alpha2)/(2 * gamma)*(1+lambda+gamma)/sqrt(Q)
  1/(2 * gamma) * 
    (alpha2 * (gamma - 1 + sqrt(Q) - 
                 lambda * (1+lambda+gamma)/sqrt(Q)) + 
       gamma * (1+(1+lambda+gamma)/sqrt(Q)))
}

alpha2=rexp(1);gamma=rexp(1);lambda=rexp(1)
c(risk(alpha2,gamma,lambda),risk2(alpha2,gamma,lambda))

## converges to 1 as alpha2 -> Inf
lazy_term <- function(alpha2, gamma, const) {
  lambda <- (alpha2 + 1)*const
  Q <- (gamma - lambda - 1)^2 + 4 * gamma * lambda
  sqrt(Q) - lambda * (1+lambda+gamma)/sqrt(Q)
}

lazy_ratio <- function(alpha2, gamma, const) 
  risk(alpha2, gamma, const * (alpha2+1))/opt_risk(alpha2, gamma)

# "lazy term -> (gamma + 1)"
# 0.1 %>% {
#   c(lazy_term(1, gamma, .),
#   lazy_term(10, gamma, .),
#   lazy_term(100, gamma, .),
#   lazy_term(1000, gamma, .),
#   lazy_term(1e4, gamma, .))  
# }

# "lazy ratio -> (gamma)/(gamma - 1)"
# gamma <- rexp(1)+1
#  5.1 %>% {
#    c(lazy_ratio(1, gamma, .),
#    lazy_ratio(10, gamma, .),
#   lazy_ratio(100, gamma, .),
#   lazy_ratio(1000, gamma, .),
#   lazy_ratio(1e4, gamma, .))  }
# (gamma)/(gamma-1)

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

# "Check lim{alpha2} risk*/alpha2"
# library(magrittr)
# gamma <- 10
# alpha2 * 1000 %>% {c(opt_risk(., gamma)/., (gamma-1)/gamma)}

alpha2 <- rexp(1); gamma <- rexp(1)
n <- 1e3; p <- floor(gamma * n)
risk(alpha2, gamma, gamma/alpha2)
avger(ridge_error_random_I, 20, 7, TRUE, alpha2, gamma, n * gamma/alpha2, n = n, naive = TRUE)
avger(ridge_error_random_I, 20, 7, TRUE, alpha2, gamma, n * gamma/alpha2, n = n)

risk(alpha2, gamma, 2 * gamma/alpha2)
avger(ridge_error_random_I, 20, 7, TRUE, alpha2, gamma, 2*n * gamma/alpha2, n = n, naive = TRUE)
avger(ridge_error_random_I, 20, 7, TRUE, alpha2, gamma, 2*n * gamma/alpha2, n = n)


risk(alpha2, gamma, 3 * gamma/alpha2)
avger(ridge_error_random_I, 20, 7, TRUE, alpha2, gamma, 3*n * gamma/alpha2, n = n, naive = TRUE)
avger(ridge_error_random_I, 20, 7, TRUE, alpha2, gamma, 3*n * gamma/alpha2, n = n)


ridge_error_random_I(alpha2, gamma, gamma/alpha2, naive = TRUE)
source("decision_theory//asymptotic_general.R")
res <- risk_formula(1, 1, alpha2, gamma)
plot(res$lambdas[res$lambdas > 0], res$risk[res$lambdas > 0], type ="l", xlim = c(0, 10))
abline(v = gamma/alpha2)
ind <- which(res$lambdas > .5 & res$lambdas < 2)[1]
(lambda <- res$lambdas[ind])
res$risk[ind]
risk(alpha2, gamma, lambda)
avger(ridge_error_random_I, 20, 7, TRUE, alpha2, gamma, lambda, naive = TRUE)

opt_risk(alpha2, gamma)
opt_risk(alpha2, gamma, TRUE)
