library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
library(glmnet, warn.conflicts = FALSE)
source('utils/zattach.R')
source('transfer/source.R')
source('eb_ident/eigenprism.R')
source('eb_ident/bayes_reg.R')
source('eb_ident/source.R')

#set.seed(1)
mcc <- 7
hyperpars <- list(n=30, pY= 100, pX=100 , W_X= 2, s_e= 1, s_b = 0.05, 
                  W_e= 2, rho_t = 0.9, L= 100, n_te= 100)
pars <- do.call(gen_params, hyperpars)
truth <- do.call(gen_data, pars)
true_e <- modifyList(truth, list(Sigma_t = eye(hyperpars$n)))
obs <- do.call(obs_data, truth)
#zattach(truth)

## Estimate parameters
pre_Bayes <- do.call2(predictive_Bayes, truth, mc.cores = mcc)
bayes_cl <- do.call(post_probs, c(obs, pre_Bayes))$cl
(bayes_err <- sum(bayes_cl != truth$i_chosen))

## bayes with wrong Sigma_t
pre_Bayes_we <- do.call2(predictive_Bayes, true_e, mc.cores = mcc)
bayes_we_cl <- do.call(post_probs, c(obs, pre_Bayes_we))$cl
(bayes_we_err <- sum(bayes_we_cl != truth$i_chosen))

## Oracle ML (with Bayes estimates)
pre_MLEb <- do.call(pre_mle, c(obs, pre_Bayes))
MLEb_cl <- do.call(post_probs, c(obs, pre_MLEb))$cl
(MLEb_err <- sum(MLEb_cl != truth$i_chosen))

p_CV <- do.call2(params_CV1, obs, filtr = FALSE, mc.cores = mcc, rule = "lambda.min")
pre_CV <- do.call(pre_mle, c(obs, p_CV))
CV_cl <- do.call(post_probs, c(obs, pre_CV))$cl
(CV_err <- sum(CV_cl != truth$i_chosen))

p_CV0 <- do.call2(params_CV1, obs, filtr = FALSE, mc.cores = mcc, rule = "lambda.1se",
                  est.t = FALSE, shrinkage = 0.5)
pre_EB <- do.call2(predictive_EP, c(obs, p_CV0), mc.cores = mcc)
EB_cl <- do.call(post_probs, c(obs, pre_EB))$cl
(EB_err <- sum(EB_cl != truth$i_chosen))

## Oracle EB
pre_oEB <- do.call2(predictive_EP, true_e, mc.cores = mcc)
oEB_cl <- do.call(post_probs, c(obs, pre_oEB))$cl
(oEB_err <- sum(oEB_cl != truth$i_chosen))

pre_oEB <- do.call2(predictive_EP, c(obs, p_CV0), Sigma_t = truth$Sigma_t, mc.cores = mcc)
oEB_cl <- do.call(post_probs, c(obs, pre_oEB))$cl
(oEB_err <- sum(oEB_cl != truth$i_chosen))
pre_CV <- do.call(pre_mle, c(obs, pre_oEB))
CV_cl <- do.call(post_probs, c(obs, pre_CV))$cl
(CV_err <- sum(CV_cl != truth$i_chosen))
