## DEVELOPMENT OF lineId package
## temporary file to be moved later

library(devtools)
library(roxygen2)
setwd('~/github/fmri/')
roxygenize('lineId')
load_all('lineId')
build('lineId')
install('lineId')
library(lineId)

setwd('lineId/vignettes/')
knitr::knit2html('identification.Rmd')

help(obs_data_filter)

pars <- gen_params()
dat <- do.call(gen_data, pars)
obs <- do.call(obs_data_filter, dat)
obsp <- modifyList(obs, list(i_chosen = dat$i_chosen))

filter_method = no_filter
filter_params = list()
forward_method = fit_ridge_kernel
forward_params = list()
Sigma_e_method = residual_offdiag
Sigma_e_params = list()
Sigma_t_method = assume_iid
Sigma_t_params = list()
backward_method = pre_mle
backward_params = list()
scoring_method = topk_score
scoring_params = list()

mc.cores = 3
res <- list(X = X, Y = Y, X_te = X_te, y_star = y_star,
            i_chosen = i_chosen, mc.cores = mc.cores)
res$filt <- do.call(filter_method, modifyList(res, filter_params))
res$B <- do.call(forward_method, modifyList(res, forward_params))
res$Sigma_e <- do.call(Sigma_e_method, modifyList(res, Sigma_e_params))
res$Sigma_t <- do.call(Sigma_t_method, modifyList(res, Sigma_t_params))
res$pre_moments <- do.call(backward_method, modifyList(res, backward_params))
res$plikes <- do.call(post_likes, res)
res$score <- do.call(scoring_method, modifyList(res, scoring_params))
res

res <- do.call(identification_pipeline1, obsp)

help(tkron_d_kron)


a1 <- 5; a2 <- 3; b1 <- 3; b2 <- 4
A <- pracma::randn(a1, a2)
B <- pracma::randn(b1, b2)
cc <- rnorm(a1 * b1)
ans <- tkron_d_kron(A, B, cc)
ansn <- t(A %x% B) %*% diag(cc) %*% (A %x% B)
f2(ans, ansn)

a1 <- 5; a2 <- 3; b1 <- 3; b2 <- 4
A <- pracma::randn(a1, a2)
B <- pracma::randn(b1, b2)
cc <- rnorm(a2 * b2)
ans <- kron_v(A, B, cc)
ansn <- (A %x% B) %*% cc
f2(ans, ansn)

n <- 100; pX <- 20; pY <- 30
X <- randn(n, pX)
B <- randn(pX, pY)
Sigma_b <- eye(pY)
Sigma_e <- cor(randn(3 * pY, pY))
Sigma_t <- cor(randn(3 * n, n))
E <- sqrtm(Sigma_t) %*% randn(n, pY) %*% sqrtm(Sigma_e)
Y <- X %*% B + E
res <- post_moments(X, Y, Sigma_e, Sigma_b, Sigma_t, TRUE, TRUE)
f2(res$Mu, B)
f2(solve(t(X) %*% X, t(X) %*% Y), B)

####
##  Testing Bayesian methods
####







