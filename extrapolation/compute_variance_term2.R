## compute the variance of prediction
## using rank statistics

## Define B_{k, d} as the matrix
## B_{ij} = E[U^(j-1)] given U ~ Beta(k-1-(i-1), i)
## which is equal to Beta(k-1-(i+1)+(j-1),i)/Beta(k-1-(i-1),i)
## We want to compute
##  norm^2(B* (B_{k, 2k, d}^T B_{k, 2k, d})^{-1} B_{k, 2k, d}^T)/sqrt(k)

# first row of coefficient_matrix2
coefficient_matrix1 <- function(k, d) t((k - 1)/(k + (0:d) - 1))

coefficient_matrix2 <- function(k, d) {
  ans <- matrix(0, k - 1, d + 1)
  ii <- row(ans); jj <- col(ans)
  ans <- matrix(beta(k-1-(ii-1)+(jj-1), ii)/beta(k-1-(ii-1), ii), nrow = nrow(ans))
  ans
}

vfunc2 <- function(k, K, d, rcoeffs = FALSE) {
  B <- coefficient_matrix2(k, d)
  coeffs <- coefficient_matrix1(K, d) %*% pinv(t(B) %*% B) %*% t(B)
  if (rcoeffs) return(coeffs)
  #sum(coeffs^2)
  sum(abs(coeffs))^2
}

check_unbiased2 <- function(k, K, d) {
  B <- coefficient_matrix2(k, d)
  BK <- coefficient_matrix1(K, d)
  #print(pinv(t(B) %*% B) %*% t(B) %*% B)
  rbind(BK, BK %*% pinv(t(B) %*% B) %*% t(B) %*% B)
}

## plot of v-functions

# plot((3:15), sapply(3:15, vfunc2, K = 40, d = 1), 
#      xlab = expression(k[max]), ylab = "v", ylim = c(0, 60),
#      main = "K = 40, d = 1 to 3", type = "o", lwd = 2)
# lines((4:15), sapply(4:15, vfunc2, K = 40, d = 2), col = "green", lwd = 2, type = "o")
# lines((6:15), sapply(6:15, vfunc2, K = 40, d = 3), col = "red", lwd = 2, type = "o")
# legend(10, 40, legend = c("d = 1", "d = 2", "d = 3"), 
#        col = c("black", "green", "red"),
#        lwd = 2, pch = "o")
# 
# plot((10:30), sapply(10:30, vfunc2, K = 100, d = 2), 
#      xlab = expression(k[max]), ylab = "v", ylim = c(0, 60),
#      main = "K = 100, d = 2 to 5", type = "l", lwd = 2)
# lines((10:30), sapply(10:30, vfunc2, K = 100, d = 3), col = "green", lwd = 2)
# lines((10:30), sapply(10:30, vfunc2, K = 100, d = 4), col = "red", lwd = 2)
# lines((10:30), sapply(10:30, vfunc2, K = 100, d = 5), col = "blue", lwd = 2)
# legend(45, 50, legend = c("d = 2", "d = 3", "d = 4", "d = 5"), 
#        col = c("black", "green", "red", "blue"),
#        lwd = 2)
