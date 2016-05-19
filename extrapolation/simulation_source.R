## Generate data
library(pracma)
gen.data <- function(p, sigma, K, r1, r2) {
  mus <- randn(K, p) # cluster centroids
  Z1 <- rep(1:K, each = r1)
  Z2 <- rep(1:K, each = r2)
  Ytr <- mus[Z1, ] + sigma * randn(K * r1, p) # final training data
  Yte <- mus[Z2, ] + sigma * randn(K * r2, p) # final test data
  list(Ytr = Ytr, Yte = Yte, Z1 = Z1, Z2 = Z2)  
}

## Train models

library(glmnet)
library(MASS)
library(kknn)
library(nnet)
epsilon_nn <- 0.1 # epsilon for epsilon-NN
## get log probs for first k classes 
pred_submodel <- function(Ytr, Yte, Z1, Z2, k) {
  Ytr_sub <- Ytr[Z1 %in% 1:k, ] # subset training data
  Yte_sub <- Yte[Z2 %in% 1:k, ] # subset test data
  Z1_sub <- Z1[Z1 %in% 1:k]
  Z2_sub <- Z2[Z2 %in% 1:k]
  ## gmm
  mu_hat <- t(sapply(1:k, function(i) {
    colMeans(Ytr_sub[Z1_sub == i, ])
  }))
  dist_gmm <- pdist2(mu_hat, Yte_sub)
  pred_gmm <- -t(dist_gmm)
  ## QDA
  res_qda <- qda(x = Ytr_sub, grouping = Z1_sub)
  pred_qda <- predict(res_qda, Yte_sub)$posterior
  ## multinomial logistic
  res_glmnet <- glmnet(Ytr_sub, Z1_sub, family = "multinomial")
  pred_glmnet <- predict(res_glmnet, Yte_sub, s = 0)[, , 1]
  ## eps-NN
  df_train <- data.frame(Z = as.factor(Z1_sub), Y = Ytr_sub)
  df_test <- data.frame(Z = as.factor(Z2_sub * 0 + 1), Y = Yte_sub)
  res_enn <- kknn::kknn(Z ~ ., train = df_train, test = df_test, 
                        k = floor(epsilon_nn * length(Z1_sub)))
  pred_enn <- res_enn$prob
  ## nnet
  res_nnet <- nnet(Z ~ ., data = df_train, size = 10, trace = FALSE)
  pred_nnet <- predict(res_nnet, df_test)
  list(pred_gmm = pred_gmm, pred_qda = pred_qda, pred_glmnet = pred_glmnet,
       pred_enn = pred_enn, pred_nnet = pred_nnet)
}

## Get Vij

get_vij <- function(pred, Z2) {
  rankconv <- t(apply(pred, 1, function(v) rank(v, ties.method = "random")))
  Z2_sub <- Z2[1:nrow(pred)]
  Vs <- rankconv[cbind(1:nrow(pred), Z2_sub)]
  Vs
}

## Exponential extrapolation

library(nnls)
expmix <- function(ws, as, xs) {
  as.numeric(ws %*% exp(t(t(as)) %*%  t(xs)))
}
expbasis <- function(as, xs) {
  t(exp(t(t(as)) %*%  t(xs)))
}
fit_expmix <- function(as, xs, y) {
  X <- expbasis(as, xs)
  res <- nnls(X, y)
  sol <- res$x
  sol[sol < 1e-10] <- 0
  fit_a <- as[sol > 0]
  fit_w <- sol[sol > 0]
  ff <- function(xs) expmix(fit_w, fit_a, xs)
  list(a = fit_a, w = fit_w, f = ff)
}
binmom <- function(succ, tot, k) {
  choose(succ, k)/choose(tot, k)
}
expmix_binmom <- function(Vs, k, K, plot = FALSE) {
  momks <- sapply(1:k, function(x) mean(binmom(Vs, k, x - 1)))
  res <- fit_expmix(-seq(0, 5, 0.01), 1:k, momks)
  if (plot) {
    plot(1:max(K), res$f(1:max(K)), type = "l", ylim = c(0, 1))
    points(1:k, momks)
  }
  res$f(K)
}

## Pseudolikelihood

library(nloptr)
us = seq(0, 1, 0.02) # change the discretization level

fit_pm_models <- function(Ys, k, us = seq(0, 1, 0.02), gu_init = rep(1/length(us), length(us)),
                          K = k) {
  Ys <- as.numeric(Ys)
  (ws <- sapply(0:k, function(i) sum(Ys == i)))
  (momK <- mean(binmom(Ys, k, k)))
  usk <- us^k
  binprobs <- matrix(0, k + 1, length(us))
  for (i in 1:length(us)) binprobs[, i] <- dbinom(0:k, k, us[i])
  of_gu <- function(gu) {
    ft <- binprobs %*% gu
    -sum(ws * log(ft))
  }
  gof <- function(gu) {
    ft <- binprobs %*% gu
    -as.numeric(t(binprobs) %*% (ws/ft))
  }
  ## MPLE unconstrained
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = function(gu) sum(gu) - 1, 
                eval_jac_g_ineq = function(gu) 0 * gu + 1, 
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,
                            check_derivatives = FALSE,maxeval = 1e4))
  (t2u <- proc.time() - t1)
  # print(res)
  gu_mple <- res$solution
  # disp_solution(res$solution); title("uncon")
  
  
  ## MPLE moment constraint
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = function(gu) sum(gu) - 1, 
                eval_jac_g_ineq = function(gu) 0 * gu + 1, 
                eval_g_eq = function(gu) sum(gu * usk) - momK,
                eval_jac_g_eq = function(gu) usk,
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,
                            check_derivatives = FALSE,maxeval = 1e4))
  (t2mom <- proc.time() - t1)
  # print(res)
  gu_mom <- res$solution
  
  ## MPLE monotonic constraint
  ll <- length(us)
  mat <- matrix(0, ll - 1, ll)
  cmat <- (row(mat) == col(mat)) - (row(mat) == (col(mat) - 1))
  eval_g_ineq_mono = function(gu) c(sum(gu) - 1, gu[-ll] - gu[-1])
  eval_jac_g_ineq_mono = function(gu) rbind(0 * gu + 1, cmat)
  
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = eval_g_ineq_mono, 
                eval_jac_g_ineq = eval_jac_g_ineq_mono, 
                # eval_g_eq = function(gu) sum(gu * usk) - momK,
                # eval_jac_g_eq = function(gu) usk,
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,
                            check_derivatives = FALSE,maxeval = 1e4))
  (t2mono <- proc.time() - t1)
  # print(res)
  gu_mono <- res$solution

  ## MPLE 2 constraint
  ll <- length(us)
  mat <- matrix(0, ll - 1, ll)
  cmat <- (row(mat) == col(mat)) - (row(mat) == (col(mat) - 1))
  eval_g_ineq_mono = function(gu) c(sum(gu) - 1, gu[-ll] - gu[-1])
  eval_jac_g_ineq_mono = function(gu) rbind(0 * gu + 1, cmat)
  
  t1 <- proc.time()
  res <- nloptr(gu_init, of_gu, eval_grad_f = gof, 
                lb = 0 * us, ub = 0 * us + 1, 
                eval_g_ineq = eval_g_ineq_mono, 
                eval_jac_g_ineq = eval_jac_g_ineq_mono, 
                eval_g_eq = function(gu) sum(gu * usk) - momK,
                eval_jac_g_eq = function(gu) usk,
                opts = list(algorithm = "NLOPT_LD_SLSQP",
                            xtol_rel = 1.0e-8,
                            print_level = 0,check_derivatives = FALSE,
                            maxeval = 1e4))
  (t2c2 <- proc.time() - t1)
  # print(res)
  gu_mm <- res$solution
  
  list(gu_mple = gu_mple, gu_mom = gu_mom, gu_mono = gu_mono, gu_mm = gu_mm,
       mple_est = sum(gu_mple * us^(K  - 1)),
       mom_est = sum(gu_mom * us^(K  - 1)),
       mono_est = sum(gu_mono * us^(K  - 1)),
       mm_est = sum(gu_mm * us^(K  - 1)))
}

## High-dimensional asymptotics

meanexp <- function (v) 
{
    vm <- max(v)
    mean(exp(v - vm)) * exp(vm)
}

piK <- function (mus, K, mc.reps = 10000) 
{
    samp <- qnorm(((1:mc.reps) - 0.5)/mc.reps)
    sampmat <- repmat(t(samp), length(mus), 1) - mus
    temp <- log(1 - pnorm(sampmat))
    1 - apply((K - 1) * temp, 1, meanexp)
}

inv_piK <- function (p, K, upper = 20, res = 1000) 
{
    if (p == 0) 
        return(Inf)
    xs <- seq(0, upper, length.out = res + 1)
    ps <- piK(xs, K)
    xs[order(abs(ps - p))[1]]
}

extrapolate_hd <- function(acck, k, K) {
  1-piK(inv_piK(1-acck, k), K)
}

## Build table

extrapolation_methods <- c("acc_sub", "acc_full", "exp", 
                           "pmle", "pmle_mom", "pmle_mono", "pmle_mm",
                           "hd")
classifiers <- c("gmm", "qda", "glmnet", "enn", "nnet")

build_extrapolation_table <- function(synth_data, k, K) {
  
  tab <- matrix(NA, length(classifiers), length(extrapolation_methods), 
                dimnames = list(classifiers, extrapolation_methods))
  preds_sub <- pred_submodel(synth_data$Ytr, synth_data$Yte, 
                           synth_data$Z1, synth_data$Z2, k)
  preds_full <- pred_submodel(synth_data$Ytr, synth_data$Yte, 
                            synth_data$Z1, synth_data$Z2, K)
  V_subs <- lapply(preds_sub, get_vij, Z2 =synth_data$Z2)
  acc_sub <- sapply(preds_sub, function(v) mean(get_vij(v, synth_data$Z2) == k))
  acc_full <- sapply(preds_full, function(v) mean(get_vij(v, synth_data$Z2) == K))
  tab[, "acc_full"] <- acc_full
  tab[, "acc_sub"] <- acc_sub
  tab[, "exp"] <- sapply(V_subs, expmix_binmom, k = k, K = K)
  pmle_fits <- lapply(V_subs, fit_pm_models, k = k, us = us, K = K)
  for (i in 1:length(pmle_fits)) {
    tab[i, c("pmle", "pmle_mom", "pmle_mono", "pmle_mm")] <-
      unlist(pmle_fits[[i]][5:8])
  }
  tab[, "hd"] <- sapply(acc_sub, extrapolate_hd, k = k, K = K)
  tab
}

