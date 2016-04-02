####
##  Predictions with bootstrap
##  get `P' distribution
####

source("Yuval/ident_setup_minus.R")
ms <- c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250)
nm <- length(ms)

## create training split
load("Yuval/v1_data.RData")
standards = apply(fit_feat,2,sd)
constF = which(standards==0)     

for (nvox in c(20, 40, 60, 80, 100)) {

fs <- paste0("Yuval/replicates/pred_rules_boots_", 100:109, ".RData")
Pemps <- matrix(0, length(fs), 250)
for (iii in 1:length(fs)) {
  load(fs[iii])
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
  usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]
  
  ##trainpred = getPreds(c_trainF,ind_struct=vox_prediction_rules,voxind=usevox)
  testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
  validpred = getPreds(c_validF,ind_struct=vox_prediction_rules,voxind=usevox)
  
  newCondVar =genCondVar(1:10409,usevox)
  
  # make variance symmetric
  newCondVarS = (newCondVar + t(newCondVar))/2
  ridgeinv = solve(newCondVarS + diag(length(usevox))*3)
  
  scores <- matchResults(testpred, testY, 'MSE',invCov = ridgeinv, rankit = FALSE)
  Pemp <- 1 - sapply(1:250, function(i) sum(scores[i, -i] < scores[i, i])/249)
  Pemps[iii, ] <- Pemp
}
fname <- paste0('Yuval/replicates/Pemp_res_nvox', nvox, '.rds')
saveRDS(Pemps, file = fname)
}

