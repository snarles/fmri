## approximating larger bump nonnegatively with smaller bumps that are regularly spaced

library(nnls)


shift <- function(v, places) {
  if (length(places) > 1) return(sapply(places, function(p) shift(v, p)))
  if (places == 0) return(v)
  if (places > 0)
    return(c(v[(places+1):length(v)], v[1:places]))
  if (places < 0)
    return(c(v[(length(v) + places + 1):length(v)], v[1:(length(v) + places)]))
}

delt <- 0.01
master.grid <- seq(-10, 10, by = delt)

subgrid <- seq(-5, 5, by = 0.25)

basis_mat <- shift(dnorm(master.grid, sd = 0.5), subgrid/delt)

matplot(master.grid, basis_mat, type = "l")

residual_calc <- function(mu, sd = 1, plot = FALSE) {
  signal <- dnorm(master.grid, mean = mu, sd)
  nnls_res <- nnls(basis_mat, signal)
  if(plot) {
    plot(master.grid, signal, type = "l")
    lines(master.grid, nnls_res$fitted, col = "red")
  }
  max(abs(nnls_res$residuals))
}

pts <- sort(runif(100, min = -1.2, max = 1.2))
resids <- sapply(pts, residual_calc)

plot(pts, resids, type = "o")
hist(resids)
