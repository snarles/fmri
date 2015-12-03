source("topological/boot_btb_source.R")

get_p_value <- function(XA, XB, G0, G1, GA, GB, SigmaA, SigmaB,
                        int.frac = 0, rep.out = 100,
                        rep.in = 100, mcc = 0) {
  BA0 <- G0 %*% GA
  BB0 <- ((1-int.frac) * G0 + int.frac * G1) %*% GB
  YA <- XA %*% BA0 + mvrnorm(nA, rep(0, q), SigmaA)
  YB <- XB %*% BB0 + mvrnorm(nB, rep(0, q), SigmaA)
  vals <- double_boot(boot_stat, XA, YA, XB, YB, rep.out = rep.out,
                      rep.in = rep.in, mcc = mcc)
  veczero_pvalue(vals)
}

####
## p-values
####

## parameter setup

rep.out <- 500
rep.in <- 200
mcc <- 39

nA <- 30; nB <- 30; p <- 2; q <- 5
XA <- randn(nA, p)
XB <- randn(nB, p)
G0 <- randn(p, q); G1 <- randn(p, q)
GA <- svd(randn(q, q))$u
GB <- svd(randn(q, q))$u
SigmaA <- cov(randn(2 * q, q))
SigmaB <- cov(randn(2 * q, q))

## computation

int.fracs <- c(0, 0.05, 0.1, 0.2)
neach <- 100
factors <- rep(1:length(int.fracs), each = neach)
pvals <- 0* factors

for (i in 1:length(factors)) {
  pvals[i] <- get_p_value(XA, XB, G0, G1, GA, GB, SigmaA, SigmaB,
                          int.fracs[factors[i]], rep.out, rep.in, mcc)
}

## plotting results

pdf("boot_btb_plot1.pdf")

layout(matrix(1:4, 2, 2))
for (i in 1:length(int.fracs)) {
  hist(pvals[factors == i], main = paste("effect size", int.fracs[i]))
}

dev.off()

pdf("boot_btb_plot2.pdf")

layout(1)
plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1), ylab = "p-value", xlab = "quantile")
for (i in 1:length(int.fracs)) {
  lines(1:neach/neach, sort(pvals[factors==i]))
}
abline(0, 1, col = "red")
title("effect sizes 0, 0.05, 0.1, 0.2")

dev.off()
