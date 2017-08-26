## L-infinity polynomial approximation
## us = grid points, u_i
## fu = function values, f(u_i)
## d = degree

library(linprog)
library(pracma)

uniform_poly <- function(d, us, fu) {
  n <- length(us)
  #plot(us, fu, type = "o")
  
  Vmat <- repmat(t(t(us)), 1, d + 1) ^ repmat(0:d, n, 1)
  
  ## L infinity norm
  x.t <- c(t = 0)
  ## Polynomial coefficients
  x.ap <- rep(0, d + 1); names(x.ap) <- paste0("a", 0:d, "p")
  x.an <- rep(0, d + 1); names(x.an) <- paste0("a", 0:d, "n")
  ## Slack variables
  x.sp <- rep(0, n); names(x.sp) <- paste0("s", 1:n, "p")
  x.sn <- rep(0, n); names(x.sn) <- paste0("s", 1:n, "n")
  xvec <- c(x.t, x.ap, x.an, x.sp, x.sn)
  
  ## Polynomial equation
  Amat.poly <- cbind(0, Vmat, -Vmat, eye(n), -eye(n))
  colnames(Amat.poly) <- names(xvec)
  rownames(Amat.poly) <- paste0("poly", 1:n)
  #View(Amat.poly)
  const.poly <- rep("==", nrow(Amat.poly))
  b.poly <- fu
  names(b.poly) <- rownames(Amat.poly)
  names(const.poly) <- rownames(Amat.poly)
  
  ## Slack constraint
  Amat.slack <- cbind(-1, zeros(2 * n, 2 * (d + 1)), eye(2 * n))
  colnames(Amat.slack) <- names(xvec)
  rownames(Amat.slack) <- names(c(x.sp, x.sn))
  #View(Amat.slack)
  const.slack <- rep("<=", nrow(Amat.slack))
  b.slack <- rep(0, nrow(Amat.slack))
  names(b.slack) <- rownames(Amat.slack)
  names(const.slack) <- rownames(Amat.slack)
  
  ## combine the equality and inequality constraints
  Amat <- rbind(Amat.poly, Amat.slack)
  bvec <- c(b.poly, b.slack)
  const.dir <- c(const.poly, const.slack)
  #View(data.frame(const.dir, bvec, Amat))
  
  ## Call linprog!!
  #help(solveLP)
  res <- linprog::solveLP(c(1, rep(0, length(xvec) - 1)),
                          bvec, Amat, const.dir = const.dir, lpSolve = TRUE)
  res$opt
  ss <- res$solution
  names(ss) <- names(xvec)
  coeffs <- ss[1 + (1:(d+1))] - ss[d + 2 + (1:(d+1))]
  fu_tilde <- Vmat %*% coeffs
  #c(max(abs(fu_tilde - fu)), res$opt)
  #plot(us, fu, type = "o"); lines(us, fu_tilde, col = "red")
  list(coeffs = coeffs, error = res$opt, 
       fu_tilde = fu_tilde, us = us, fu = fu)
}

# d <- 7
# us <- seq(-1, 1, 0.01)
# fu <- sin(5 * us)
# res <- uniform_poly(d, us, fu)
# with(res, {
#   plot(us, fu, type = "o"); lines(us, fu_tilde, col = "red")
# })
# res$error
