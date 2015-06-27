gen_params <- function(n=30, pY= 60, pX= 70, W_X= 2, s_e= 10, W_e= 2, L= 100, n_te= 10, ...) {
  Sigma_X <- 1/(W_X * pX) * randn(W_X * pX, pX) %>% { t(.) %*% . }
  Sigma_e <- s_e/(W_e * pY) * randn(W_e * pY, pY) %>% { t(.) %*% . }
  s0s <- 1* rnorm(pY)^2 + 1
  Sigma_b <- diag(s0s)
  dSigma_B <- rep(s0s, each = pX)
  Bvec <- sqrt(dSigma_B) * rnorm(pY * pX)
  B0 <- matrix(Bvec, pX, pY)
  X <- mvrnorm(n, mu = rep(0, pX), Sigma = Sigma_X)
  X_te <- mvrnorm(L, mu = rep(0, pX), Sigma = Sigma_X)
  list(X = X, X_te = X_te, B = B0, Sigma_b = Sigma_b, Sigma_e = Sigma_e, n_te = n_te)  
}

gen_data <- function(X, X_te, B, Sigma_b, Sigma_e, n_te, ...) {
  n <- dim(X)[1]
  L <- dim(X_te)[1]
  pX <- dim(X)[1]
  pY <- dim(B)[2]
  Y <- X %*% B + mvrnorm(n, mu = rep(0, pY), Sigma = Sigma_e)
  i_chosen <- sample(L, n_te, TRUE)
  y_star <- X_te[i_chosen, , drop = FALSE] %*% B +  mvrnorm(n_te, rep(0, pY), Sigma_e)
  list(X = X, Y = Y, X_te = X_te, i_chosen = i_chosen, y_star = y_star, B = B,
       Sigma_b = Sigma_b, Sigma_e = Sigma_e)
}

obs_data <- function(X, X_te, Y, y_star, ...) list(X = X, X_te = X_te, Y = Y, y_star = y_star)

predictive_Bayes <- function(X, Y, X_te, Sigma_e, Sigma_b, ...) {
  pX <- dim(X)[2]; pY <- dim(Y)[2]
  B <- post_moments(X, Y, Sigma_e, Sigma_b, computeCov = FALSE)
  res <- post_predictive(X, Y, X_te, Sigma_e, Sigma_b)
  list(pre_moments = res, filt = rep(TRUE, pY), B = B, Sigma_e = Sigma_e)
}

post_probs <- function(X_te, y_star, pre_moments, filt, ...) {
  y_filt <- y_star[, filt]
  L <- dim(X_te)[1]
  n_te <- dim(y_filt)[1]
  pY <- dim(y_filt)[2]
  pprobs <- matrix(0, n_te, L)
  for (i in 1:L) {
    Mu <- pre_moments[[i]]$Mu
    Cov <- pre_moments[[i]]$Cov
    ld <- -log(det(Cov))
    resid <- t(t(y_filt) - Mu)
    ss <- solve(Cov, t(resid))
    ips <- rowSums(resid * t(ss))
    pprobs[, i] <-  -ld - ips
  }
  cl <- apply(pprobs, 1, function(v) order(-v)[1])
  list(pprobs = pprobs, cl = cl)
}

params_CV1 <- function(X, Y, ...) {
  n <- dim(X)[1]
  pX <- dim(X)[2]
  pY <- dim(Y)[2]
  B <- matrix(0, pX, pY) # placeholder for estimate
  resids <- matrix(0, n, pY)
  lambdas <- numeric(pY) # lambda for each column of B
  for (i in 1:pY) {
    res <- cv.glmnet(X, Y[, i], alpha = 0, intercept = FALSE, standardize = FALSE,
                     grouped = FALSE)
    pre <- predict(res, newx = X, s = res$lambda.min)
    resids[, i] <- pre - Y[, i]
    B[, i] <- coef(res, newx = X, s = res$lambda.min)[-1]
    lambdas[i] <- (Norm(B[, i])^2/pX)  
  }
  filt <- (lambdas > 1e-5)
  lambdas <- lambdas[filt]
  B <- B[, filt, drop = FALSE]
  resids <- resids[, filt, drop = FALSE]
  Sigma_e <- 0.5 * cov(resids) + 0.5 * diag(diag(cov(resids)))
  list(pre_moments = NULL, filt = filt, B = B, Sigma_e = Sigma_e, Sigma_b = diag(lambdas))
}

pre_mle <- function(X_te, B, Sigma_e, filt, ...) {
  L <- dim(X_te)[1]
  ans <- as.list(numeric(L))
  Mus <- X_te %*% B
  for (i in 1:L) ans[[i]] <- list(Mu = Mus[i, ], Cov = Sigma_e)
  list(pre_moments = ans, filt = filt, B = B, Sigma_e = Sigma_e)
}








# params_EP <- function(X, Y, Xte, Yte, avg = FALSE) {
#   q <- dim(X)[2]
#   p <- dim(Y)[2]
#   nte <- dim(Yte)[1]
#   
#   res_EP <- eigenprisms(X, Y)
#   T2 <- res_EP$T2
#   if (avg) {
#     T2 <- rep(mean(T2), length(T2))
#   }
#   filt_EP <- (T2 > 0)
#   lambdas_EP <-  T2[filt_EP]/q
#   p_EP <- sum(filt_EP)
#   B_mu_EP <- matrix(0, q, p_EP)
#   Y_EP <- Y[, filt_EP]
#   for (i in 1:p_EP) {
#     B_mu_EP[, i] <- solve(t(X) %*% X + diag(rep(lambdas_EP[i], q))) %*% t(X) %*% Y_EP[, i]    
#   }
#   Yhat_EP <- X %*% B_mu_EP
#   Sigma_e_EP <- cov(Y_EP - Yhat_EP) %>% {0.5 * . + 0.5 * diag(1e-4 + diag(.))}
#   Sigma_B_EP <- diag(rep(1/lambdas_EP, each = q))
#   EP_cov_B_vec <- solve(solve(Sigma_e_EP) %x% (t(X) %*% X) + diag(1/diag(Sigma_B_EP)))
#   list(lambdas = lambdas_EP, filt = filt_EP, B = B_mu_EP, Sigma_e = Sigma_e_EP,
#        Sigma_B = EP_cov_B_vec, res_EP = res_EP)
# }
