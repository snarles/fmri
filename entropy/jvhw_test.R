source("entropy//jvhw.R")

ss <- 30
pdist <- runif(ss)
pdist <- pdist/sum(pdist)

(h_true <- sum(-pdist * log(pdist)))

n <- 1e2
samp <- sample(1:ss, n, replace = TRUE, prob = pdist)

freqs <- table(samp)
(h_m <- h_mle(freqs))
(h_j <- h_jvhw(freqs))

c(h_true, h_m, h_j)
