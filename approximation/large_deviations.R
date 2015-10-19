source("approximation/large_deviations_source.R")

sigma2 <- 1
ds <- c(1:10, 100 * 1:30, 1e4 * 1:10)
ress <- matrix(0, length(ds), 3)
for (i in 1:length(ds)) {
  d <- ds[i]
  theta <- sqrt((1 + sigma2) * d); r <- sqrt(sigma2 * d)
  ress[i, ] <-   c(log_true_ans(d, theta, r),
  log_prox_ans(d, theta, r),
  log_prox_ans2(d, theta, r))  
}

View(cbind(ds, ress))




d=1200
theta=48.98979
r=34.64102
x=20.28261

deriv2_logFp(d, theta, r, x)


ff <- dbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
FF <- pbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
fp <- deriv_dbeta((1-ip_x(theta, r, x)) / 2, (d-1)/2, (d-1)/2)
di <- deriv_ip_x(theta, r, x)
di2 <- deriv2_ip_x(theta, r, x)

log(di)
log(ff) -log(FF) + log(-di)

x=0.1464467
a=b=599.5
deriv_dbeta(x, a, b)
