####
##  Predictions with bootstrap
####

getPreds = function(dataF,ind_struct=vox_prediction_rules,voxind=usevox){
  preds = matrix(0,nc = length(voxind),nr = nrow(dataF))
  for (i in 1:length(voxind)){
    preds[,i] = dataF %*% ind_struct[[voxind[i]]]$beta_coef
  }
  return (preds)
}

genCondVar = function(featind,voxind,ind_struct=vox_prediction_rules){
  nvox = length(voxind)
  voxcov = (t(trainY[,voxind])%*%trainY[,voxind])/(nrow(trainY))
  voxfeatcov = (t(trainY[,voxind]) %*% c_trainF[,featind])/(nrow(trainY))
  fullbeta = matrix(0,nc=length(voxind),nr=length(featind))
  for (i in 1:length(voxind)){
    fullbeta[,i] = ind_struct[[voxind[i]]]$beta_coef
  }
  voxCondCov = voxcov - voxfeatcov%*%fullbeta
  return (voxCondCov)
}

bestMatch = function(preds, dataY, score,invCov = diag(length(usevox)), voxind=usevox){
  npoints = nrow(preds)
  chosen = numeric(nrow(dataY))
  for (i in 1:nrow(dataY)){
    if (score=='MSE'){
      Yi = matrix(rep(dataY[i,voxind],npoints),nr=npoints,byrow=T)
      ResYi = Yi - preds
      chosen[i] = which.min(diag(ResYi %*% invCov %*% t(ResYi)))
    }
    else if (score=='COR'){
      chosen[i] = which.max(apply(preds,1,cor,dataY[i,voxind]))
    }
    else if (score=='COV'){
      chosen[i] = which.max(apply(preds,1,cov,dataY[i,voxind]))
    }
    if ((i%%20)==0){
      cat(',')
    }
  }
  return(sum(chosen==(1:nrow(dataY)))/nrow(dataY))
}

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
