####
##  Extrapolation of information curves
####

## verification of some theoretical calculations

library(lineId)
pk <- function(iota, ...) 1 - piK(sqrt(2 * iota), ...)



iota <- 0.5

pk(iota, 2)
pk(iota, 3)
pk(iota, 4)

us <- 1:99999/100000
gg <- dnorm(qnorm(us) - sqrt(2 * iota))/dnorm(qnorm(us))
sum(gg)/length(us)
plot(us, gg, type = "l")
gg2 <- gg/sum(gg)

sum(us * gg2)
sum(us^2 * gg2)
sum(us^3 * gg2)
