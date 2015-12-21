set.seed(0)
n <- 900
samp <- sample(1:99, n, TRUE)
for (i in 1:15) samp <- sample(samp, n, TRUE)
samp <- sort(samp)
samp <- c(samp, rep(100, 100))

table(table(samp))
39/n

# for (i in 1:(length(samp)/20)) {
#   cat(paste0(paste0(samp[((i-1)* 20+1) : (i * 20)], collapse = ";"), ";"))
# }

