## finite-sample classification models with generated data
## function returns a list with:
##  - plikes : matrix of margins
##  - i_chosen : true class labels
## true class = row number

library(pracma)
library(MASS)
library(magrittr)

simulate_qda <- function(
  k = 100,
  r_train = 20,
  r_test = 1,
  p = 2,
  wishart_par = 3, # larger -> covariances closer to identity (Inf for identity)
  snr = 3, 
  Sigma_mu = snr * eye(p),
  bayes = FALSE # Bayes oracle
) {
  ## Generate true parameters
  mus <- mvrnorm(k, rep(0, p), Sigma_mu)
  if (wishart_par < Inf) {
    Sigmas <- lapply(1:k, function(ind) {
      x <- randn(wishart_par * p, p)
      (t(x) %*% x)/(wishart_par * p)
    })
    xerrs <- lapply(1:k, function(ind) {
      mvrnorm(r_train + r_test, rep(0, p), Sigmas[[ind]])
    })
  } else {
    Sigmas <- rep(list(eye(p)), k)
    xerrs <- lapply(1:k, function(ind) randn(r_train + r_test, p))
  }
  ## Estimate the parameters from data
  if (bayes) {
    mu_hats <- mus
    Sigma_hats <- Sigmas
  } else {
    mu_errs <- lapply(xerrs, function(x) colMeans(x[(1:r_train), ])) %>%
      do.call(what = rbind)
    mu_hats <- mus + mu_errs
    Sigma_hats <- lapply(xerrs, function(x) cov(x[(1:r_train), ]))
  }
  iSigma_hats <- lapply(Sigma_hats, solve)
  Sigma_dets <- sapply(Sigma_hats, det)
  iSigma_stack <- lapply(iSigma_hats, as.numeric) %>% do.call(what = rbind)
  xtest <- lapply(1:k, function(ind) {
    repmat(mus[ind, ], r_test, 1) + xerrs[[ind]][-(1:r_train), ]
  }) %>%
    do.call(what = rbind)
  i_chosen <- rep(1:k, each = r_test)
  ## Compute the matrix of margins, dimension is
  ##  (est_class k * true_class k * repeat r_test) x p 
  diffs <- repmat(xtest, k, 1) - (mu_hats %x% ones(nrow(xtest), 1))
  ## form quadratic expansion, dimension is
  ##  (est_class k * true_class k * repeat r_test) x p^2
  diffs2 <- repmat(diffs, 1, p) * (diffs %x% ones(1, p))
  ssquares <- rowSums(diffs2 * (iSigma_stack %x% ones(nrow(xtest), 1))) %>%
    matrix(nrow = nrow(xtest))
  ## check correctness:
  i_ind <- sample(k, 1);j_ind <- sample(k, 1);r_ind <- sample(r_test, 1)
  rbind(diffs[((j_ind - 1) * k + (i_ind - 1))*r_test + r_ind, ],
    xtest[i_ind, ] - mu_hats[j_ind, ])
  c(ssquares[(i_ind-1)*r_test + r_ind, j_ind],
    t(xtest[i_ind, ] - mu_hats[j_ind, ]) %*% iSigma_hats[[j_ind]] %*%
      (xtest[i_ind, ] - mu_hats[j_ind, ]))
  plikes <- -1/2 * (ssquares + repmat(log(Sigma_dets), nrow(xtest), 1))
  ## return answer
  list(plikes = plikes, i_chosen = i_chosen)
}

# res <- simulate_qda(snr = 20, wishart_par = 1, bayes = TRUE)
# # res <- simulate_qda(snr = 20, wishart_par = Inf, bayes = TRUE)
# do.call(lineId::topk_scores, res)
# source("approximation/gaussian_lc_source.R")
# 1 - mc(20 * eye(2), 100)
