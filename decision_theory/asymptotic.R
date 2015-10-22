####
##  Asymptotic ridge
####

gamma <- 2.5
mI <- function(z) {
  -((z + gamma - 1) + sqrt((z + gamma-1)^2 - 4 * gamma * z))/
    (2 * gamma * z)
}
mId <- function(z) {
  -mI(z)/z + 1/(2 * gamma * z) *
    (-1 - (z- gamma - 1)/sqrt((z + gamma - 1)^2 - 4 * gamma*z))
}

mI_num <- function(z) ((z + gamma - 1) + sqrt((z + gamma-1)^2 - 4 * gamma * z))
mI_num_d <- function(z) 1 + (z-gamma-1)/sqrt((z + gamma-1)^2-4*gamma*z)

mI_num(-1e-10)
c(numDeriv::grad(mI_num, -2),mI_num_d(-2))
c(mI_num_d(-1e-3),-2/(gamma - 1))

c(numDeriv::grad(mI, -1),mId(-1))

mI(-1e-10)
-1/(gamma)/(gamma - 1)
