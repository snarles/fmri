source("topological/boot_btb_source.R")

get_p_value <- function(XA, XB, G0, G1, GA, GB, SigmaA, SigmaB,
                        int.frac = 0, rep.out = 100, mcc = 0) {
  BA0 <- G0 %*% GA
  BB0 <- ((1-int.frac) * G0 + int.frac * G1) %*% GB
  YA <- XA %*% BA0 + mvrnorm(nA, rep(0, q), SigmaA)
  YB <- XB %*% BB0 + mvrnorm(nB, rep(0, q), SigmaA)
  vals <- double_boot(boot_stat, XA, YA, XB, YB, rep.out = rep.out, mcc = mcc)
  veczero_pvalue(vals)
}

####
## p-values
####

rep.out <- 10
mcc <- 7

nA <- 40; nB <- 50; p <- 2; q <- 10
XA <- randn(nA, p)
XB <- randn(nB, p)
G0 <- randn(p, q); G1 <- randn(p, q)
GA <- svd(randn(q, q))$u
GB <- svd(randn(q, q))$u
SigmaA <- cov(randn(2 * q, q))
SigmaB <- cov(randn(2 * q, q))

int.fracs <- rep(c(0, 0.05, 0.1), each = 10)
pvals <- 0* int.fracs

for (i in 1:length(int.fracs)) {
  pvals[i] <- get_p_value(XA, XB, G0, G1, GA, GB, SigmaA, SigmaB,
                          int.fracs[i], rep.out, mcc)
}

pdf("boot_btb_plot1.pdf")
layout(1)
boxplot(pvals ~ int.fracs, data.frame(int.fracs, pvals), main = "pvalues")
dev.off()