####
##  Predictions with bootstrap
####


## Prepare data

library(lineId)

matchResults = function(preds, dataY, score,invCov = diag(length(usevox)), voxind=usevox, rankit = TRUE){
  npoints = nrow(preds)
  res <- list()
  for (i in 1:nrow(dataY)){
    if (score=='MSE'){
      Yi = matrix(rep(dataY[i,voxind],npoints),nr=npoints,byrow=T)
      ResYi = Yi - preds
      res[[i]] = diag(ResYi %*% invCov %*% t(ResYi))
    }
    else if (score=='COR'){
      res[[i]] = -apply(preds,1,cor,dataY[i,voxind])
    }
    else if (score=='COV'){
      res[[i]] = -apply(preds,1,cov,dataY[i,voxind])
    }
    if ((i%%20)==0){
      cat(',')
    }
  }
  res <- do.call(rbind, res)
  if (rankit) res <- t(apply(res, 1, rank))
  res
}

source('Yuval/reducedModels/identification_functions.R')

load('Yuval/v1_data.RData')
standards = apply(fit_feat,2,sd)
constF = which(standards==0)     

ms <- c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250)
nm <- length(ms)

## random voxel selection
vox1 = order(SNRv1_corr,decreasing = TRUE)[1:1250]
set.seed(1)
vox1 <- sample(vox1, 1250, FALSE)

for (nvox in c(100, 200, 300, 400, 500, 600)) {
##for (nvox in c(20, 40, 60, 80)) {
fs <- paste0("Yuval/reducedModels/replicates/pred_rules_sel5_boots_", 100:109, ".RData")
ihats <- matrix(0, length(fs), nm)
for (iii in 1:length(fs)) {
  load(fs[iii], verbose = TRUE)
  fit_trainS <- attr(vox_prediction_rules_sel5, "train_samp")
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
  
  max_pyr_level = 5
  
  # generates the feature attributes vector.
  # Orientation -> featAtt[1,] (1 vertical, 5 horizontal)  
  # Pyramid level -> featAtt[2,]
  # Vertical location -> featAtt[3,]  (vert/horz might be confused)
  # Horizontal location -> featAtt[4,]
  
  featAll = featAtt()
  featNoConst = featAll[,-constF]
  select_feat = which(featNoConst[2,] <= max_pyr_level)
  
  sel_trainF = c_trainF[,select_feat]
  sel_testF = c_testF[,select_feat]
  sel_validF = c_validF[,select_feat]
  validY = v1$val
  
  # choose best nvox voxels
  ## usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]
  usevox = vox1[1:nvox]
  
    
  ##trainpred = getPreds(c_trainF,ind_struct=vox_prediction_rules,voxind=usevox)
  sel_trainpred = getPreds(sel_trainF,ind_struct=vox_prediction_rules_sel5,usevox)
  sel_testpred = getPreds(sel_testF,ind_struct=vox_prediction_rules_sel5,usevox)
  sel_validpred = getPreds(sel_validF,ind_struct=vox_prediction_rules_sel5,usevox)
  #imagepred = getPreds(c_aimagesF,ind_struct=voxel_predicts,voxind=usevox)
  
  newCondVar_sel =genCondVar(dataF = sel_trainF,dataY = trainY,1:ncol(sel_trainF),usevox,ind_struct = vox_prediction_rules_sel5)
  
  # make variance symmetric
  newCondVarS_sel = (newCondVar_sel + t(newCondVar_sel))/2
  ridgeinv_sel = solve(newCondVarS_sel + diag(length(usevox))*3)
  
  scores <- matchResults(sel_testpred, testY, 'MSE',invCov = ridgeinv_sel, rankit = FALSE)
  ihat <- numeric(nm)
  for (ind in 1:nm) {
    m <- ms[ind]
    (p <- resample_misclassification(-scores, 1:250, m, replace = FALSE))
    (IH <- Ihat_LI(p, m, 20))
    ihat[ind] <- IH
  }
  ihats[iii, ] <- ihat
}
fname <- paste0('Yuval/reducedModels/replicates/res_voxseed01_nvox', nvox, '.rds')
saveRDS(ihats, file = fname)
}

