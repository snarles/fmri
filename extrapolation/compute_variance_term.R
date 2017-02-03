## compute the variance of prediction

## Define B_{k1, k2, d} as the matrix
## B_{ij} = (k1 + i - 2)/(k1 + i + j - 3) for j = 1,..,(d+1) and i = 1,...,k2-k1+1
## We want to compute
##  norm^2(B_{K, K, d} (B_{k, 2k, d}^T B_{k, 2k, d})^{-1} B_{k, 2k, d}^T)/sqrt(k)

coefficient_matrix <- function(k1, k2, d) {
  ans <- matrix(0, k2 - k1 + 1, d + 1)
  ii <- row(ans); jj <- col(ans)
  ans <- matrix((k1 + ii - 2)/(k1 + ii + jj - 3), nrow = nrow(ans))
  ans
}

vfunc <- function(k, K, d, rcoeffs = FALSE) {
  B <- coefficient_matrix(k, 2*k, d)
  coeffs <- coefficient_matrix(K, K, d) %*% pinv(t(B) %*% B) %*% t(B)
  if (rcoeffs) return(coeffs)
  sum(coeffs^2)
}

check_unbiased <- function(k, K, d) {
  B <- coefficient_matrix(k, 2*k, d)
  BK <- coefficient_matrix(K, K, d)
  #print(pinv(t(B) %*% B) %*% t(B) %*% B)
  rbind(BK, BK %*% pinv(t(B) %*% B) %*% t(B) %*% B)
}

## plot of v-functions

plot(2*(2:15), sapply(2:15, vfunc, K = 40, d = 1), 
     xlab = expression(k[max]), ylab = "v", ylim = c(0, 50),
     main = "K = 40, d = 1 to 3", type = "o", lwd = 2)
lines(2*(5:15), sapply(5:15, vfunc, K = 40, d = 2), col = "green", lwd = 2, type = "o")
lines(2*(7:15), sapply(7:15, vfunc, K = 40, d = 3), col = "red", lwd = 2, type = "o")
legend(20, 40, legend = c("d = 1", "d = 2", "d = 3"), 
       col = c("black", "green", "red"),
       lwd = 2, pch = "o")

plot(2*(10:30), sapply(10:30, vfunc, K = 100, d = 2), 
     xlab = expression(k[max]), ylab = "v", ylim = c(0, 60),
     main = "K = 100, d = 2 to 5", type = "l", lwd = 2)
lines(2*(10:30), sapply(10:30, vfunc, K = 100, d = 3), col = "green", lwd = 2)
lines(2*(10:30), sapply(10:30, vfunc, K = 100, d = 4), col = "red", lwd = 2)
lines(2*(10:30), sapply(10:30, vfunc, K = 100, d = 5), col = "blue", lwd = 2)
legend(45, 50, legend = c("d = 2", "d = 3", "d = 4", "d = 5"), 
       col = c("black", "green", "red", "blue"),
       lwd = 2)
