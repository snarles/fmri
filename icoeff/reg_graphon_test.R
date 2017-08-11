source("icoeff/reg_graphon.R")

p <- 5
max.pow = 5
P <- gen_rwg(p = p, self.loop = FALSE)
S <- gen_rwg(p = p, self.loop = FALSE)

rbind(pspec(P, max.pow = max.pow), 
      pspec(P %*% S, max.pow = max.pow), 
      pspec(S %*% P, max.pow = max.pow))
