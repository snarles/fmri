####
##  Check the glm coefficients
####

source('Yuval/ident_setup.R')

###
## Run the identification with pre-computed rules
###
nvox = 1250# up to 1250
usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]
testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
newCondVar =genCondVar(1:10409,usevox)
# make variance symmetric
newCondVarS = (newCondVar + t(newCondVar))/2
ridgeinv = solve(newCondVarS + diag(length(usevox))*3)

bestMatch(testpred,testY,'MSE',invCov = ridgeinv)

####
##  Check the glmnet coefficients
####

library(glmnet)
vox.ind <- 100
rule <- vox_prediction_rules[[vox.ind]]
feat.inds <- which(rule$beta_coef != 0)
fit <- glmnet(c_trainF, trainY[, vox.ind], alpha = rule$alpha)
coef <- coef(fit, s = rule$lambda)
plot(coef[-1], rule$beta_coef)

fit <- glmnet(c_trainF[, feat.inds], trainY[, vox.ind], alpha = rule$alpha)
coef <- coef(fit, s = rule$lambda)
plot(coef[-1], rule$beta_coef[feat.inds])
