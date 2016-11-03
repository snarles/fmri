library(MASS)
library(pracma)
library(glmnet)
library(lineId)

load('Yuval/v1_data.RData')
load('Yuval/pred_rules.RData')

standards = apply(fit_feat,2,sd)
constF = which(standards==0)     

# Choose training set
load('Yuval/fit_trainS.RData') # loading fit_trainS

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

getPreds = function(dataF,ind_struct=vox_prediction_rules,voxind=usevox){
  preds = matrix(0,nc = length(voxind),nr = nrow(dataF))
  for (i in 1:length(voxind)){
    preds[,i] = dataF %*% ind_struct[[voxind[i]]]$beta_coef
  }
  return (preds)
}

# choose best nvox voxels
#nvox = 1250# up to 1250
nvox <- 200
usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]


trainpred = getPreds(c_trainF,ind_struct=vox_prediction_rules,voxind=usevox)
testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
validpred = getPreds(c_validF,ind_struct=vox_prediction_rules,voxind=usevox)
#imagepred = getPreds(c_aimagesF,ind_struct=voxel_predicts,voxind=usevox)

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
newCondVar =genCondVar(1:10409,usevox)

# make variance symmetric
newCondVarS = (newCondVar + t(newCondVar))/2
ridgeinv = solve(newCondVarS + diag(length(usevox))*3)

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


scores <- matchResults(testpred, testY, 'MSE',invCov = ridgeinv, rankit = FALSE)

binom_lcb <- function(acc, nt, alpha = 0.05) {
  qbinom(alpha, nt, acc, lower.tail = TRUE)/nt
}

source("mi_vs_aba/var_resample.R")

alpha <- 0.05
(zalpha <- -qnorm(alpha, lower.tail = TRUE))

(p <- resample_misclassification(-scores, 1:250, 250, replace = FALSE))
acc <- 1-p
aba_to_mi_lower(250, acc)
Ihat_LI(1-acc, 250)

(acc_l <- binom_lcb(acc, 250, alpha) - 1/sqrt(4 * 250))
aba_to_mi_lower(250, acc_l)

k <- 5

ks <- 2:250

ils <- c(); res2 <- c()
for (k in ks) {
  K <- 250
  (p <- resample_misclassification(-scores, 1:250, k, replace = FALSE))
  acc <- 1-p
  # ve <- var_est(-scores, k)
  # sqrt(ve)
  acc_l <- acc - zalpha * sqrt((1-acc)*(acc)/K) - 1/sqrt(4 * K)
  ils <- c(ils, aba_to_mi_lower(k, acc_l))
  #res2 <- c(res2, Ihat_LI(p, k))
}

plot(ks, ils, type = "o")




