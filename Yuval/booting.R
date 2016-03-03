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
v1$val <- v1$val[, 1:1250]

rownames(v1$resp) <- casenames
rownames(fit_feat) <- casenames
colnames(v1$resp) <- voxnames
colnames(v1$val) <- voxnames


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

makedata2 <- function(fit_trainS, sel_inds, voxind) {
  trainF = fit_feat[fit_trainS,-constF][, sel_inds]
  trainY = v1$resp[fit_trainS,][, voxind]
  testF = fit_feat[-fit_trainS,-constF][, sel_inds]
  testY = v1$resp[-fit_trainS,][, voxind]
  validF = val_feat[,-constF][, sel_inds]
  validY = v1$val[, voxind]
  
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

get_boot <- function(i) {
  set.seed(i)
  train_inds <- sort(sample(1750, 1500))
  zattach(makedata(train_inds))
  vox_prediction_rules = lapply(voxinds, run_glmnet)
  names(vox_prediction_rules) <- voxinds
  # trainpred = getPreds(c_trainF,ind_struct=vox_prediction_rules,voxind=voxinds)
  testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=voxinds)
  list(rules = vox_prediction_rules, testpred = testpred)
}

####
##  Actual bootstraps for a few voxels
####

voxinds <- order(SNRv1_corr,decreasing = TRUE)[1:9] 
  ##vox179 vox474 vox396 vox826 vox711 vox691 vox493 vox658 vox811
voxinds <- colnames(v1$resp)[voxinds]

# testpreds <- vector("list", 1750); 
# rules <- list()
# names(testpreds) <- casenames
# for(i in 1:1750) testpreds[[i]] <- matrix(0, 0, length(voxinds))
# library(parallel)
# boot.reps <- 40
# t1 <- proc.time()
# res <- mclapply(40 + (1:boot.reps), get_boot, mc.cores = 40)
# proc.time() - t1
# for (ii in 1:boot.reps) {
#   rules[[ii]] <- res[[ii]]$rules
#   testpred <- res[[ii]]$testpred
#   for (i in rownames(testpred)) {
#     testpreds[[i]] <- rbind(testpreds[[i]], testpred[i, ])
#   }
# }
# 
# save(testpreds, rules, file = "Yuval/temp_booting2.RData")

load("Yuval/aggregate_booting.RData")
rules <- rules_all; rm(rules_all)
testpreds <- testpreds_all; rm(testpreds_all)
gc()

####
##  Try to do the bootstraps
####

(vox <- voxinds[4])
ruleind <- 2
# (caseind <- casenames[3])
# rules[[ruleind]][[vox]]$lambda
# rules[[ruleind]][[vox]]$alpha
# sum(rules[[ruleind]][[vox]]$beta_coef != 0)
feats <- which(rules[[ruleind]][[vox]]$beta_coef != 0)
mc.reps <- 120
preds <- matrix(NA, 1750, mc.reps, dimnames = list(casenames))
for (i in 1:mc.reps) {
  train_inds <- sort(sample(1750, 1500))
  zattach(makedata2(train_inds, feats, vox))
  pred <- c_testF %*% solve(t(c_trainF) %*% c_trainF) %*% t(c_trainF) %*% trainY
  preds[rownames(pred), i] <- pred
}

#hist(preds[caseind, ])
#hist(testpreds[[caseind]][, vox])

mus <- lapply(testpreds, function(v) mean(v[, vox]))
sds <- lapply(testpreds, function(v) sd(v[, vox]))
musP <- apply(preds, 1, mean, na.rm = TRUE)
sdsP <- apply(preds, 1, sd, na.rm = TRUE)

##View(cbind(mus, musP, sds, sdsP))
layout(t(1:2))
plot(mus, musP); abline(0, 1, col = "red")
plot(sds, sdsP); abline(0, 1, col = "red")