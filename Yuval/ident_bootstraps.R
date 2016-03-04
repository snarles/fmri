####
##  Charles' stuff using Yuval's stuff
##  Bootstrap the accuracy numbers
####

source('Yuval/ident_setup.R')
standards = apply(fit_feat,2,sd)
constF = which(standards==0)   
ttF <-  fit_feat[, -constF]
c_ttF = scale(ttF,center = TRUE, scale = TRUE)
ttY = v1$resp[, 1:1250]

## Use the same condvar to save time (!?)
# choose best nvox voxels
nvox = 500# up to 1250
usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]
testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
validpred = getPreds(c_validF,ind_struct=vox_prediction_rules,voxind=usevox)
newCondVar =genCondVar(1:10409,usevox)
# make variance symmetric
newCondVarS = (newCondVar + t(newCondVar))/2
ridgeinv = solve(newCondVarS + diag(length(usevox))*3)

## make random CV splits
mc.reps <- 10
nte <- 250
train_splits <- t(apply(t(1:100), 2, function(v) {
  sample(c(rep(TRUE, 1750-nte), rep(FALSE, nte)), 1750, FALSE)  
}))
times <- c()
results <- matrix(NA, mc.reps, 1)
for (split.ind in 1:mc.reps) {
  t1 <- proc.time()
  tr.inds <- train_splits[split.ind, ]
  tpred.new <- matrix(0, nte, nvox)
  for (vox.no in 1:nvox) {
    vox.ind <- usevox[vox.no]
    rule <- vox_prediction_rules[[vox.ind]]
    feat.inds <- which(rule$beta_coef != 0)
    if (length(feat.inds) > 0) {
      c_ttSub <- c_ttF[, feat.inds]
      teY <- ttY[!tr.inds, vox.ind]
      trY <- ttY[tr.inds, vox.ind]
      teX <- c_ttSub[!tr.inds, ]
      trX <- c_ttSub[tr.inds, ]
      Yh <- teX %*% ginv(trX) %*% trY
      tpred.new[, vox.no] <- Yh
    }
  }
  score1 <- bestMatch(tpred.new, ttY[!tr.inds, ],'MSE',invCov = ridgeinv)
  results[split.ind, ] <- c(score1)
  (times[split.ind] <- (proc.time() - t1)["elapsed"])
}


bestMatch(validpred,validY[1:120,],'MSE',invCov = ridgeinv)
bestMatch(testpred,testY[1:120,],'MSE',invCov = ridgeinv)
bestMatch(testpred,testY[1:250,],'MSE',invCov = ridgeinv)