####
##  Does the min of ncx converge to a distribution?
####

## log 1+x
log_1_plus <- function(x) {
  if (abs(x) > 1e-3) return(log(1 + x))
  sum(-1 * (-x)^(1:20)/(1:20))
}

## med_loc breaks down at k > 20
med_loc <- function(k, ncp, L) qchisq(1 - (.5)^(1/L), k, ncp)

med_loc2 <- function(k, ncp, L, lx = 0, ux = 100, nits = 100) {
  const <- -log(2)
  for (i in 1:nits) {
    xc <- (ux + lx)/2
    val <- L * log_1_plus(-pchisq(xc, k, ncp))
    if (val > const) {
      lx <- xc
    } else {
      ux <- xc
    }
    #print(xc)
  }
  xc
}


get_L <- function(x, k, ncp) {
  pp <- pchisq(x, k, ncp)
  log(1/2)/log_1_plus(-pp)
}

#x <- 0.1; k <- 1; ncp <- 1.1
#x=60;k=600;ncp=660
get_log_L <- function(x, k, ncp) {
  pp <- pchisq(x, k, ncp, log.p = TRUE)
  if (pp > -10) {
    pp <- log(-log_1_plus(-exp(pp)))
  }
  ans <- log(log(2)) - pp
  if (ans == Inf) cat(paste0("x=", x, ";k=", k, ";ncp=", ncp))
  ans
}

mincdf <- function(x, k, z, L) 1 - exp(L * log_1_plus(-pchisq(x, k, k + z * sqrt(k))))
minpdf <- function(x, k, z, L) {
  pp <- pchisq(x, k, k + z * sqrt(k))
  L * exp(L * log_1_plus(-pp))/(1 - pp)
}

# get_L(1, 100, -2.5)
# 
# 
# 
# k <- 30
# L <- get_L(1, k, 0)
# med_loc(k, 0, L)
# med_loc2(k, 0, L)
# 
# 
# Ls <- matrix(0, 100, 10)
# for (k in 1:10) {
#   Ls[, k] <- sapply(1:100, function(i) get_L(k, i, 0))  
# }
# matplot(log(Ls), type = "l")
# 
# 
# 
# 
# 
# 
# 
# mincdf(4, 100, Ls[100, 5])
# mincdf(5, 100, Ls[100, 5])
# mincdf(6, 100, Ls[100, 5])
# 
# minpdf(5, 80, Ls[80, 5])
# minpdf(5, 100, Ls[100, 5])
# 
# 
# ## what happens with differnt y^2?
# 
# 
# x_med <- 5
# Xs <- matrix(0, 91, 10)
# lalas <- -5:4/2
# for (k in 1:10) {
#   Xs[, k] <- sapply(10:100, function(i) med_loc2(i, lalas[k], Ls[i, x_med]))  
# }
# matplot(10:100, Xs, type = "l")

####
##  Fix sigma2, increase d, fix error prob
####

fixed_error_L <- function(d, sigma2) {
  ey2 <- (1 + sigma2) * d
  alpha <- sigma2/(1 + sigma2)
  ed1 <- alpha * (alpha * ey2 + d)
  get_log_L(ed1, d, ey2)
}

sigma2 <- 3
ds <- 1:1000
logLs <- sapply(ds, function(d) fixed_error_L(d, sigma2))
plot(logLs, type = "l")




