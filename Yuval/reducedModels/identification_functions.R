featAtt = function(){
  # featAttVec: row1 - orientation, row2 - pyr-level, row3 - locationx, row4 - lovationy
  pyr  = c(1,8,4*8,16*8,64*8,256*8,1024*8)
  lens = c(0,1,2,   4,   8,   16,    32,  64)
  featAttVec = matrix(0,nr=4,nc=sum(pyr))
  featAttVec[1,] = (0:(sum(pyr)-1)) %% 8
  cumsumpyr = cumsum(pyr)
  featAttVec[2,1] = 1;
  featAttVec[3:4,1] = 0;
  for (i in 1:(length(cumsumpyr)-1)) {
    k = 0
    for (j in seq((cumsumpyr[i]+1),cumsumpyr[i+1],8)) {
      featAttVec[2,j:(j+7)] = i+1
      featAttVec[3,j:(j+7)] = floor(k/(lens[i+1])) + 1
      featAttVec[4,j:(j+7)] = (k%%(lens[i+1])) + 1
      k = k + 1
    }
  }
  return(featAttVec)
}


require(glmnet)
run_glmnet = function(voxnum, dataF, dataY, alpha_val=0.8, nfolds = 6) {
  res = cv.glmnet(dataF, dataY[,voxnum],nfolds = 6,alpha = alpha_val)  
  best_lam = res$lambda.min
  best_ind = which(res$lambda==best_lam)
  pred_rule = res$glmnet.fit$beta[,best_ind]
  pred_score = res$cvm[best_ind]
  return(list(lambda = best_lam, alpha = alpha_val,beta_coef= pred_rule, score=pred_score))
}

getPreds = function(dataF,ind_struct=vox_prediction_rules,voxind=usevox){
  preds = matrix(0,nc = length(voxind),nr = nrow(dataF))
  for (i in 1:length(voxind)){
    preds[,i] = dataF %*% ind_struct[[voxind[i]]]$beta_coef
  }
  return (preds)
}


genCondVar = function(dataF, dataY,featind,voxind,ind_struct){
  nvox = length(voxind)
  voxcov = (t(dataY[,voxind])%*%dataY[,voxind])/(nrow(dataY))
  voxfeatcov = (t(dataY[,voxind]) %*% dataF[,featind])/(nrow(dataY))
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

