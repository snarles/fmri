####
##  Spherical cap area
####

library(hypergeo)

log_vsph <- function(d, r=1)
  log(2/d) + (d/2)*log(pi) - lgamma(d/2) + d*log(r)

mc.s <- 1e6
p <- 3
## [-1, 1]^p hypercube
pts <- 2*matrix(runif(mc.s * p) - .5, mc.s, p)
nms <- rowSums(pts^2)
ball <- pts[nms < 1, ]

## check volume
c(mean(nms < 1) * 2^p, exp(log_vsph(p)))

## cap probability
emp_prob <- function(u) mean(ball[, 1] > (1- 2*u))
integral_prob <- function(u, res = 1e3) {
  xs <- seq(-1, 1, length.out=res)
  vslice <- (sqrt(1-xs^2))^(p-1)
  sum(vslice[xs > (1-2*u)])/sum(vslice)
}
theory_prob <- function(u) {
  ff <- function(x) x * hypergeo(1/2, (1-p)/2, 3/2, x^2)
  Re((ff(1) - ff(1-2*u))/(2*ff(1)))
}

c(0.5, emp_prob(0.5), integral_prob(0.5), theory_prob(0.5))

runif(1) %>% {c(emp_prob(.), integral_prob(.), theory_prob(.))}

