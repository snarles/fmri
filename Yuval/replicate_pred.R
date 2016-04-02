####
##  Predictions with bootstrap
####

source("Yuval/ident_setup_minus.R")
ms <- c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250)
nm <- length(ms)

## create training split
load("Yuval/v1_data.RData")
standards = apply(fit_feat,2,sd)
constF = which(standards==0)     
load("Yuval/replicates/pred_rules_boots_101.RData")
fit_trainS <- attr(vox_prediction_rules, "train_samp")
trainF = fit_feat[fit_trainS,-constF]
trainY = v1$resp[fit_trainS,]
testF = fit_feat[-fit_trainS,-constF]
testY = v1$resp[-fit_trainS,]
validF = val_feat[,-constF]
validY = v1$val
c_trainF = scale(trainF,center = TRUE, scale = TRUE)
feat_mean_vec = attr(c_trainF,'scaled:center')
feat_scale_vec = attr(c_trainF,'scaled:center')
c_testF = scale(testF,center = feat_mean_vec, scale = feat_scale_vec)
c_validF = scale(validF,center = feat_mean_vec, scale = feat_scale_vec)


# choose best nvox voxels
nvox = 500# up to 1250
usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]

##trainpred = getPreds(c_trainF,ind_struct=vox_prediction_rules,voxind=usevox)
testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
validpred = getPreds(c_validF,ind_struct=vox_prediction_rules,voxind=usevox)

newCondVar =genCondVar(1:10409,usevox)

# make variance symmetric
newCondVarS = (newCondVar + t(newCondVar))/2
ridgeinv = solve(newCondVarS + diag(length(usevox))*3)

bestMatch(validpred,validY[1:120,],'MSE',invCov = ridgeinv)
bestMatch(testpred,testY[1:120,],'MSE',invCov = ridgeinv)
bestMatch(testpred,testY[1:250,],'MSE',invCov = ridgeinv)


scores <- matchResults(testpred, testY, 'MSE',invCov = ridgeinv, rankit = FALSE)
ihat <- numeric(nm)
for (ind in 1:nm) {
  m <- ms[ind]
  (p <- resample_misclassification(-scores, 1:250, m, replace = FALSE))
  (p_prime <- rbinom(1, 250, p)/250)
  (IH <- Ihat_LI(p_prime, m, 20))
  ihats[ind] <- IH
}
