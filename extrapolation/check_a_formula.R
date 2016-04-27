#####
##  Checking out model: monotonic g model
#####

k <- 7
dgu <- runif(k + 1)
gu <- cumsum(dgu)
(gu <- gu/sum(gu))
(dgu <- c(gu[1], gu[-1] - gu[-(k + 1)]))
# rbind(cumsum(dgu), gu)
binprobs <- matrix(0, k + 1, k + 1)
for (i in 0:k) binprobs[, i + 1] <- dbinom(0:k, k, i/k)
# matplot(binprobs, type = "l")
hy <- binprobs %*% gu
vprobs <- matrix(0, k + 1, k + 1)
for (i in 0:k) vprobs[, i+1] <- binprobs %*% c(rep(0, i), rep(1, k + 1 - i))
hy2 <- vprobs %*% dgu
cbind(hy, hy2)
