#' Exponentiate then take the mean
#' 
#' @param v Values v_i
#' @return Mean of exp(v_i)
#' @export
#' @examples 
#' meanexp(c(-1, 3, 4))
#' mean(exp(c(-1, 3, 4)))
meanexp <- function(v) {
  vm <- max(v)
  mean(exp(v - vm)) * exp(vm)
}

#' Function which computes a Gaussian probability
#' 
#' Computes the probability that a normal with mean mu is smaller than the max of K-1 other iid Gaussian
#' @param mus Mean mu, can be vectorized.
#' @param K Number of independent gaussians.
#' @param mc.reps Numerical resolution (larger for more precision.)
#' @return A probability
#' @export
#' @examples 
#' piK(3, 10)
piK <- function(mus, K, mc.reps = 1e4) {
  samp <- qnorm(((1:mc.reps) - 0.5)/mc.reps)
  sampmat <- repmat(t(samp), length(mus), 1) - mus# one row per mu  
  temp <- log(1 - pnorm(sampmat))
  1 - apply((K-1) * temp, 1, meanexp)
}

#' Function which inverts piK
#' 
#' \code{piK} computes the probability that a normal with mean mu is smaller than the max of K-1 other iid Gaussian.
#' This the the inverse function.
#' @param p Probability
#' @param K Number of independent gaussians.
#' @param upper A guess on the upper bound of mu.
#' @param res Numerical resolution (larger for more precision.)
#' @return The parameter mu such that \code{piK(mu, K) = p}.
#' @export
#' @examples 
#' p <- piK(3, 10)
#' inv_piK(p, 10, 5)
inv_piK <- function(p, K, upper = 10, res = 1e3) {
  if (p ==0) return(Inf)
  xs <- seq(0, upper, length.out = res + 1)
  ps <- piK(xs, K)
  xs[order(abs(ps- p))[1]]
}

#' Low-information estimator of Mutual Information
#' 
#' Given a misclassification rate and number of classes, estimate the mutual information between X and Y.
#' @param p Misclassification probability.
#' @param K Number of classes.
#' @param upper A guess on the upper bound of the mutual information (should be exponential in K.)
#' @param res Numerical resolution (larger for more precision.)
#' @return The estimated mutual information.
#' @export
#' @examples 
#' Ihat_LI(0.5, 5, 100)
Ihat_LI <- function(p, K, upper = 10, res = 1e3) {
  1/2 * inv_piK(p, K, sqrt(2 * upper), res)^2
}

