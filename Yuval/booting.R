####
##  Bootstraps without fitting the full elastic net
##  by Charles
####

library(lineId)
load('Yuval/v1_data.RData')
standards = apply(fit_feat,2,sd)
constF = which(standards==0)    

casenames <- paste0("s", 1:1750)
voxnames <- paste0("vox", 1:1250)
v1$resp <- v1$resp[, 1:1250]
rownames(v1$resp) <- casenames
rownames(fit_feat) <- casenames
colnames(v1$resp) <- voxnames

run_glmnet = function(voxnum, alpha_val=0.8, nfolds = 6) {
  require(glmnet)
  res = cv.glmnet(c_trainF, trainY[,voxnum],nfolds = 6,alpha = alpha_val)  
  best_lam = res$lambda.min
  best_ind = which(res$lambda==best_lam)
  pred_rule = res$glmnet.fit$beta[,best_ind]
  pred_score = res$cvm[best_ind]
  return(list(lambda = best_lam, alpha = alpha_val,beta_coef= pred_rule, score=pred_score))
}

makedata <- function(fit_trainS) {
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
  list(c_trainF = c_trainF, trainY = trainY, c_testF = c_testF, c_validF = c_validF )
}

getPreds = function(dataF,ind_struct=vox_prediction_rules,voxind=usevox){
  preds = matrix(0,nc = length(voxind),nr = nrow(dataF))
  for (i in 1:length(voxind)){
    preds[,i] = dataF %*% ind_struct[[voxind[i]]]$beta_coef
  }
  rownames(preds) <- rownames(dataF)
  colnames(preds) <- voxind
  return (preds)
}

####
##  Actual bootstraps for a few voxels
####

voxinds <- order(SNRv1_corr,decreasing = TRUE)[1:9]
voxinds <- colnames(v1$resp)[voxinds]
rules <- list()
testpreds <- vector("list", 1750); 
names(testpreds) <- casenames
for(i in 1:1750) testpreds[[i]] <- matrix(0, 0, length(voxinds))
library(parallel)
boot.reps <- 20
t1 <- proc.time()
for (i in 1:boot.reps) {
  set.seed(i)
  train_inds <- sort(sample(1750, 1500))
  zattach(makedata(train_inds))
  vox_prediction_rules = mclapply(voxinds, run_glmnet, mc.cores = 3)
  names(vox_prediction_rules) <- voxinds
  rules[[i]] <- vox_prediction_rules
  ##trainpred = getPreds(c_trainF,ind_struct=vox_prediction_rules,voxind=voxinds)
  testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=voxinds)
  for (i in rownames(testpred)) {
    testpreds[[i]] <- rbind(testpreds[[i]], testpred[i, ])
  }
}
proc.time() - t1

save(testpreds, rules, file = "Yuval/temp_booting.RData")
