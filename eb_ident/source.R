mle_cl <- function(X, Y0, i_chosen, pars) {
  Y <- Y0[, pars$filt]
  L <- dim(X)[1]
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  mu_hat <- X %*% pars$B
  a <- isqrtm(pars$Sigma_e)
  mu_hat2 <- mu_hat %*% a
  Y2 <- Y %*% a
  pprobs <- (-pdist2(mu_hat2, Y2)) %>% {. - apply(., 2, min)} %>% exp %>% {./rowSums(.)}
  cl <- apply(pprobs, 1, function(v) order(-v)[1])
  (mle_err <- sum(cl != i_chosen))
  list(cl = cl, err = mle_err, pprobs = pprobs)
}

post_probs <- function(X, Y0, i_chosen, pars) {
  Y <- Y0[, pars$filt]
  L <- dim(X)[1]
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  mus <- list(L)
  covs <- list(L)
  pprobs <- matrix(0, n, L)
  for (i in 1:L) {
    mus[[i]] <- X[i, , drop = FALSE] %*% pars$B
    covs[[i]] <- (diag(rep(1, p)) %x% t(X[i, ])) %*% 
      pars$Sigma_B %*% (diag(rep(1, p)) %x% t(t(X[i, ]))) + pars$Sigma_e
    for (j in 1:n) {
      pprobs[j, i] <- -log(det(covs[[i]])) - 
        (Y[j, , drop = FALSE] - mus[[i]]) %*% solve(covs[[i]]) %*% 
        t(Y[j, , drop = FALSE] - mus[[i]])
    }  
  }
  cl <- apply(pprobs, 1, function(v) order(-v)[1])
  err <- sum(cl != i_chosen)
  list(pprobs = pprobs, cl = cl, err = err)
}

bagged_cl <- function(X, Y, i_chosen, ps) {
  kk <- length(ps)
  L <- dim(X)[1]
  nte <- dim(Y)[1]
  cl_votes <- matrix(0, nte, L)
  for (i in 1:kk) {
    res <- mle_cl(X, Y, i_chosen, ps[[i]])
    cl_votes %<>% add(res$pprobs)
  }
  cl <- apply(cl_votes, 1, function(v) order(-v)[1])
  err <- sum(cl != i_chosen)
  list(cl_votes = cl_votes, cl = cl, err = err)
}

params_CV <- function(X, Y, Xte, Yte) {
  q <- dim(X)[2]
  p <- dim(Y)[2]
  nte <- dim(Yte)[1]
  B_mu_EB <- matrix(0, q, p) # placeholder for estimate
  resids <- matrix(0, nte, p)
  lambdas <- numeric(p) # lambda for each column of B
  for (i in 1:p) {
    res <- cv.glmnet(X, Y[, i], alpha = 0, intercept = FALSE, standardize = FALSE,
                     grouped = FALSE)
    pre <- predict(res, newx = Xte)
    resids[, i] <- pre - Yte[, i]
    B_mu_EB[, i] <- as.numeric(coefficients(res))[-1]  
    lambdas[i] <- (sum(as.numeric(coefficients(res))[-1]^2)/q)  
  }
  lambdas <- pmax(lambdas, 1e-5)
  Sigma_E_hat <- 0.5 * cov(resids) + 0.5 * diag(diag(cov(resids)))
  Inv_Sigma_B <- diag(rep(1/lambdas, each = q))
  EB_cov_B_vec <- solve(solve(Sigma_E_hat) %x% (t(X) %*% X) + Inv_Sigma_B)  
  list(lambdas = lambdas, B = B_mu_EB, Sigma_e = Sigma_E_hat,  Sigma_B = EB_cov_B_vec,
       filt = rep(TRUE, p))
}

params_EP <- function(X, Y, Xte, Yte, avg = FALSE) {
  q <- dim(X)[2]
  p <- dim(Y)[2]
  nte <- dim(Yte)[1]
  
  res_EP <- eigenprisms(X, Y)
  T2 <- res_EP$T2
  if (avg) {
    T2 <- rep(mean(T2), length(T2))
  }
  filt_EP <- (T2 > 0)
  lambdas_EP <-  T2[filt_EP]/q
  p_EP <- sum(filt_EP)
  B_mu_EP <- matrix(0, q, p_EP)
  Y_EP <- Y[, filt_EP]
  for (i in 1:p_EP) {
    B_mu_EP[, i] <- solve(t(X) %*% X + diag(rep(lambdas_EP[i], q))) %*% t(X) %*% Y_EP[, i]    
  }
  Yhat_EP <- X %*% B_mu_EP
  Sigma_e_EP <- cov(Y_EP - Yhat_EP) %>% {0.5 * . + 0.5 * diag(1e-4 + diag(.))}
  Sigma_B_EP <- diag(rep(1/lambdas_EP, each = q))
  EP_cov_B_vec <- solve(solve(Sigma_e_EP) %x% (t(X) %*% X) + diag(1/diag(Sigma_B_EP)))
  list(lambdas = lambdas_EP, filt = filt_EP, B = B_mu_EP, Sigma_e = Sigma_e_EP,
       Sigma_B = EP_cov_B_vec, res_EP = res_EP)
}

p_EP_2_sd <- function(X, Y, pars) {
  q <- dim(X)[2]
  p <- dim(Y)[2]
  p_EP <- sum(pars$filt)
  lambdas_EP <- pars$lambdas
  Sigma_e_EP <- pars$Sigma_e
  Sigma_B_EP <- diag(rep(1/lambdas_EP, each = q))
  Lambda <- diag(rep(lambdas_EP, each = q))
  temp <- solve((diag(rep(1, p_EP)) %x% (t(X) %*% X)) + Lambda)
  sampling_cov_B_vec <- temp %*% (Sigma_e_EP %x% diag(rep(1, q))) %*% temp
  pars2 <- pars
  pars2$Sigma_B <- sampling_cov_B_vec
  pars2
}

params_Bayes <- function(X, Y, Sigma_e, Sigma_B, s0s) {
  q <- dim(X)[2]
  p <- dim(Y)[2]
  B_mu <- matrix(0, q, p)
  for (i in 1:p) {
    B_mu[, i] <- solve(t(X) %*% X + diag(rep(1/s0s[i], q))) %*% t(X) %*% Y[, i]    
  }
  Sigma_Bhat_vec <- solve(solve(Sigma_e) %x% (t(X) %*% X) + diag(1/diag(Sigma_B)))
  list(B = B_mu, Sigma_B = Sigma_Bhat_vec, filt = rep(TRUE, p),
       Sigma_e = Sigma_e)
}

mle_bagged <- function(X, Y) {
  p <- dim(X)[2]
  q <- dim(Y)[2]
  n <- dim(X)[1]
  k1 <- 1 + floor((0:4) * n/5)
  k2 <- floor(1:5 * n / 5)
  ps <- list(5)
  sum_e <- matrix(0, q, q)
  for (i in 1:5) {
    Xtr <- X[-(k1[i]:k2[i]), ]
    Ytr <- Y[-(k1[i]:k2[i]), ]
    Xte <- X[(k1[i]:k2[i]), ]
    Yte <- Y[(k1[i]:k2[i]), ]
    ps[[i]] <- params_CV(Xtr, Ytr, Xte, Yte)
    sum_e <- sum_e + ps[[i]]$Sigma_e/5
  }
  for (i in 1:5) ps[[i]]$Sigma_e <- sum_e
  ps
}


