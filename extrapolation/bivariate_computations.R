
rho <- 0.7
mc.reps <- 1e6

## computes Bayes error for univariate normal p(x|y); OLD VERSION
# normal_bayes_error <- function(mus, sigma2s = 1, intrange = 10, delta = 1e-3) {
#   xs <- seq(-intrange, intrange, delta)
#   pds <- lapply(mus, function(mu) dnorm(xs, mean = mu, sd = sqrt(sigma2s)))
#   maxpds <- do.call(pmax, pds)
#   ba <- 1/length(mus) * sum(maxpds) * delta
#   1 - ba
# }

## computes Bayes error for univariate normal p(x|y)
normal_bayes_error <- function(mus, sigma2s = 1) {
  smu <- sort(mus)/sqrt(sigma2s)
  md <- c(-Inf, (smu[-length(smu)] + smu[-1])/2, Inf)
  lbs <- smu - md[-1]
  ubs <- smu - md[-length(md)]
  1 - mean(pnorm(ubs) - pnorm(lbs))
}

## mus is k x N matrix
normal_bayes_error_vec <- function(mus, sigma2s = 1) {
  smu <- apply(mus, 2, sort)/sqrt(sigma2s)
  md <- rbind(-Inf, (smu[-nrow(smu), ] + smu[-1, ])/2, Inf)
  lbs <- smu - md[-1, ]
  ubs <- smu - md[-nrow(md), ]
  1 - colMeans(pnorm(ubs) - pnorm(lbs))
}





des <- list()
# des[["2"]] <- list(x = c(0, seq(0, 0.5, 0.01), 0.5), 
#                    y = c(0, rep(2, 51), 0)) ## uniform [0, 0.5]
  
for (k in 2:10) {
  t1 <- proc.time()
  bes <- normal_bayes_error_vec(rho * randn(k, mc.reps), sigma2s = 1-rho^2)
  print(proc.time() - t1)

  des[[paste(k)]] <- density(bes)
  
  plot(density(bes), main = paste(k))
  lines(density(sample(bes, mc.reps/2)), col = "red")
}

saveRDS(des, "extrapolation/bc_rho0.7.rds")
