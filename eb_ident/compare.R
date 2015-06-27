library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
library(glmnet, warn.conflicts = FALSE)
source('transfer/source.R')
source('eb_ident/eigenprism.R')
source('eb_ident/bayes_reg.R')
source('eb_ident/source.R')
library(class)

#set.seed(1)
hyperpars <- list(n=20, pY= 60, pX=60 , W_X= 2, s_e= 10, W_e= 2, L= 100, n_te= 100)
pars <- do.call(gen_params, hyperpars)
truth <- do.call(gen_data, pars)
obs <- do.call(obs_data, truth)
#zattach(truth)

## Estimate parameters
pre_Bayes <- do.call(predictive_Bayes, truth)
bayes_cl <- do.call(post_probs, c(obs, pre_Bayes))$cl
(bayes_err <- sum(bayes_cl != truth$i_chosen))

pre_MLEb <- do.call(pre_mle, c(obs, pre_Bayes))
MLEb_cl <- do.call(post_probs, c(obs, pre_MLEb))$cl
(MLEb_err <- sum(MLEb_cl != truth$i_chosen))

p_CV <- do.call(params_CV1, obs)
pre_CV <- do.call(pre_mle, c(obs, p_CV))
CV_cl <- do.call(post_probs, c(obs, pre_CV))$cl
(CV_err <- sum(CV_cl != truth$i_chosen))

