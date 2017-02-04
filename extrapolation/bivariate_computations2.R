## now with illustrations!!
set.seed(0)

plot_gaussian_densities <- function(mus, sigma2s = 1,
                                    xrange = 5, delta = 0.01, 
                                    type = c("density", "region"),
                                    col = rainbow(length(mus))) {
  type <- match.arg(type)
  plot(NA, NA, xlim = c(-xrange, xrange), ylim = c(0, .45/sqrt(sigma2s)),
       xlab = "x", ylab = "density")
  xseq <- seq(-xrange, xrange, delta)
  smu <- sort(mus)
  if (type == "density") {
    for (mu in mus) lines(xseq, dnorm(xseq, mean = mu, sd = sqrt(sigma2s)),
                          lwd = 2, col = col[smu == mu])
  } else {
    md <- c(-xrange -1, (smu[-length(smu)] + smu[-1])/2, xrange + 1)
    for (i in 1:length(smu)) {
      m1 <- md[i]
      m2 <- md[i+1]
      xs <- seq(m1, m2, delta)
      dd <- dnorm(xs, smu[i], sqrt(sigma2s))
      polygon(c(xs, m2, m1), c(dd, 0, 0), col = col[i])
    }
    for (mu in mus) lines(xseq, dnorm(xseq, mean = mu, sd = sqrt(sigma2s)),
                          lwd = 1)
  }
}

## intro plot
pdf("extrapolation/illus_example1a.pdf", width = 5, height = 3)
plot_gaussian_densities(c(-2,-0.5,2), 0.5)
dev.off()
pdf("extrapolation/illus_example1b.pdf", width = 5, height = 3)
plot_gaussian_densities(c(-2,-0.5,2), 0.5, type = "region")
dev.off()

rho <- 0.7
mc.reps <- 20

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





bess <- list()
muss <- list()

for (k in 2:10) {
  mus <- rho * randn(k, mc.reps)
  bes <- normal_bayes_error_vec(mus, sigma2s = 1-rho^2)

  bess[[paste(k)]] <- bes
  muss[[paste(k)]] <- mus
}

boxplot(bess, pch = "o", xlab = "k", ylab = "Error", ylim = c(0, 1))
stripchart(bess, vertical = TRUE, xlim = c(0, 10), add = TRUE, pch = "o", 
           col = grey(0, alpha = 0.5))


## generate the auto-box plots

k <- 3
for (ind in 1:7) {
  pdf(paste0("extrapolation/autoplots/box", k+1, "_", ind, ".pdf"),
      height = 6, width = 4)
  boxplot(bess[[k]], cex = 2, pch = "o")
  stripchart(bess[[k]][-ind], vertical = TRUE, xlim = c(0, 10), add = TRUE, pch = "o", 
             col = grey(0, alpha = 0.5), cex = 2)
  points(1, bess[[k]][ind], cex = 2, col = "red", lwd = 3)
  dev.off()
  pdf(paste0("extrapolation/autoplots/dens", k+1, "_", ind, ".pdf"),
      height = 3, width = 5)
  plot_gaussian_densities(muss[[k]][ , ind], sigma2s = 1-rho^2, type = "region")
  dev.off()
  
}
