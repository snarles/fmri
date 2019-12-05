rhos <- (1:940)/1000
plot(rhos, -1/2 * log(1 - rhos^2), type = 'l')
lines(rhos, -1/2 * log(1 - rhos^4), col='red')

rhos_s <- c(1,3,5,7,9)/10

plot(-1/2 * log(1 - rhos^2), -1/2 * log(1 - rhos^4), type = 'l', asp = 1, col = "red", xlab='I(X; Y)', ylab="I(X; X')", lwd = 2, ylim = c(-0.1, 1.0))
points(-1/2 * log(1 - rhos_s^2), -1/2 * log(1 - rhos_s^4))
abline(0, 1)
# for (rr in rhos_s) {
#   text(-1/2 * log(1 - rr^2), -1/2 * log(1 - rr^4) - 0.05, paste0(rr))
# }