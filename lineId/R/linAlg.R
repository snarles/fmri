#' Inverse square root of a matrix
#' 
#' @param m matrix
#' @export
isqrtm <- function(m) {
  res <- eigen(m)
  d <- res$values
  if (min(d) < -1e-5) warning("Negative eigenvalues in isqrtm")
  d[d < 0] <- 0
  d[d > 0] <- 1/sqrt(d[d > 0])
  v <- res$vectors
  return (v %*% diag(d) %*% t(v))
}

#' Square root of a matrix
#' 
#' @param m matrix
#' @export
sqrtm <- function(m) {
  res <- eigen(m)
  d <- res$values
  if (min(d) < -1e-5) warning("Negative eigenvalues in sqrtm")
  d[d < 0] <- 0
  d[d > 0] <- sqrt(d[d > 0])
  v <- res$vectors
  return (v %*% diag(d) %*% t(v))
}

#' Frobenius norm squared
#' 
#' @param x
#' @param y (optional)
#' @export
f2 <- function(x, y = 0) sum((x - y)^2)

