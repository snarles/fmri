####
##  Does the min of ncx converge to a distribution?
####

## log 1+x
log_1_plus <- function(x) {
  if (abs(x) > 1e-3) return(log(1 + x))
  sum(-1 * (-x)^(1:20)/(1:20))
}

med_loc <- function(k, z, L) qchisq(1 - (.5)^(1/L), k, k + z * sqrt(k))

get_L <- function(x, k, z) {
  pp <- pchisq(x, k, k + z * sqrt(k))
  log(1/2)/log_1_plus(-pp)
}

mincdf <- function(x, k, z, L) 1 - exp(L * log_1_plus(-pchisq(x, k, k + z * sqrt(k))))
minpdf <- function(x, k, z, L) {
  pp <- pchisq(x, k, k + z * sqrt(k))
  L * exp(L * log_1_plus(-pp))/(1 - pp)
}

get_L(1, 100, -2.5)


L <- get_L(1, 10, 0)
med_loc(10, 0, L)

Ls <- matrix(0, 100, 10)
for (k in 1:10) {
  Ls[, k] <- sapply(1:100, function(i) get_L(k, i, 0))  
}
matplot(log(Ls), type = "l")







mincdf(4, 100, Ls[100, 5])
mincdf(5, 100, Ls[100, 5])
mincdf(6, 100, Ls[100, 5])

minpdf(5, 80, Ls[80, 5])
minpdf(5, 100, Ls[100, 5])


## what happens with differnt y^2?

Ls <- matrix(0, 91, 10)
lalas <- -5:4/2
for (k in 1:10) {
  Ls[, k] <- sapply(10:100, function(i) get_L(5, i, lalas[k]))  
}
matplot(10:100, log(Ls), type = "l")
