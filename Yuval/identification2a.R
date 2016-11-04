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

# choose best nvox voxels
#nvox = 1250# up to 1250

nvoxs <- 1:10 * 100
scoreZ <- list()
for (nvox in nvoxs) {
  usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]
  trainpred = getPreds(c_trainF,ind_struct=vox_prediction_rules,voxind=usevox)
  testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
  validpred = getPreds(c_validF,ind_struct=vox_prediction_rules,voxind=usevox)
  #imagepred = getPreds(c_aimagesF,ind_struct=voxel_predicts,voxind=usevox)
  newCondVar =genCondVar(1:10409,usevox)
  # make variance symmetric
  newCondVarS = (newCondVar + t(newCondVar))/2
  ridgeinv = solve(newCondVarS + diag(length(usevox))*3)
  scores <- matchResults(testpred, testY, 'MSE',invCov = ridgeinv, rankit = FALSE)
  scoreZ[[paste(nvox)]] <- scores
}

save(scoreZ, file = "Yuval/scores.RData")

binom_lcb <- function(acc, nt, alpha = 0.05) {
  qbinom(alpha, nt, acc, lower.tail = TRUE)/nt
}

source("mi_vs_aba/var_resample.R")

ks <- 2:250

lbZ <- list()
iestZ <- list()

for (scores in scoreZ) {
  ils <- c(); iest <- c()
  for (k in ks) {
    K <- 250
    (p <- resample_misclassification(-scores, 1:250, k, replace = FALSE))
    acc <- 1-p
    # ve <- var_est(-scores, k)
    # sqrt(ve)
    acc_l <- acc - (zalpha + 1) * sqrt(.25/K)
    ils <- c(ils, aba_to_mi_lower(k, acc_l))
    iest <- c(iest, Ihat_LI(p, k))
  }
  
  #plot(ks, ils, type = "o")
  lbZ <- c(lbZ, list(ils))
  iestZ <- c(iestZ, list(iest))
}

names(iestZ) <- names(scoreZ)
names(lbZ) <- names(scoreZ)

save(lbZ, iestZ, file = "Yuval/lbZ.RData")

lbZ <- lbZ[1:6]
iestZ <- iestZ[1:6]

plot(ks, lbZ[[length(lbZ)]], type = "l", xlab = "k", ylab = "MI", col = "white")
for (i in c(1, 3, 5)) {
  lines(ks, lbZ[[i]], col = "black", lwd = 5)
  lines(ks, lbZ[[i]], col = rainbow(length(lbZ))[i], lwd = 3)
}
legend(150, 2.5, c("100", "300", "500"), col = "black", 
       lwd = 5)
legend(150, 2.5, c("100", "300", "500"), col = rainbow(6)[c(1,3,5)], 
       lwd = c(3,3,3), bg = NA)

plot(ks, iestZ[[length(lbZ)]], type = "l", ylim = c(0, 8))
for (i in 1:length(lbZ)) {
  lines(ks,iestZ[[i]], col = rainbow(length(lbZ))[i], lwd = 3)
}


