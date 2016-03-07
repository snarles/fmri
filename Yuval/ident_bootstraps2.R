####
##  Charles' stuff using Yuval's stuff
##  Bootstrap the accuracy numbers using the same training data
####

source('Yuval/ident_setup.R')
###
## Run the identification with pre-computed rules
###
nvox = 300# up to 1250
usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]
testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
newCondVar =genCondVar(1:10409,usevox)
# make variance symmetric
newCondVarS = (newCondVar + t(newCondVar))/2
ridgeinv = solve(newCondVarS + diag(length(usevox))*3)

bestMatch(testpred,testY,'MSE',invCov = ridgeinv)

###
## Run the identification with glmnet on pre-computed active sets
###

inds_sub <- 1:1500
inds_sub <- sample(1500, 1000, FALSE)
c_ttF <- c_trainF[inds_sub, ]
trY <- trainY[inds_sub, ]
fullbeta = matrix(0,nc=nvox,nr=ncol(c_testF))
tpred.new <- matrix(0, nrow(c_testF), nvox)
for (vox.no in 1:nvox) {
  vox.ind <- usevox[vox.no]
  rule <- vox_prediction_rules[[vox.ind]]
  feat.inds <- which(rule$beta_coef != 0)
  if (length(feat.inds) < 2) feat.inds <- c(1:5, feat.inds)
  if (length(feat.inds) > 0) {
    c_ttSub <- c_ttF[, feat.inds, drop = FALSE]
    # bt <- ginv(c_ttSub) %*% trY[, vox.ind]
    # Yh <- c_testF[, feat.inds, drop = FALSE] %*% bt
    res <- glmnet(c_ttSub, trY[, vox.ind], alpha = rule$alpha)
    bt <- coef(res, s=rule$lambda)[-1]
    Yh <- c_testF[, feat.inds, drop = FALSE] %*% bt
    fullbeta[feat.inds, vox.no] <- bt
    tpred.new[, vox.no] <- Yh
  }
}
# voxcov = (t(trY[, usevox])%*%trY[, usevox])/nrow(trY)
# voxfeatcov = (t(trY[, usevox]) %*% c_ttF)/nrow(trY)
# voxCondCov = voxcov - voxfeatcov%*%fullbeta
# newCondVarS = (newCondVar + t(newCondVar))/2
# ridgeinv = solve(newCondVarS + diag(length(usevox))*3)
(score1 <- bestMatch(tpred.new, testY, 'MSE',invCov = ridgeinv))


## make random CV splits
mc.reps <- 3
nte <- 250
train_splits <- t(apply(t(1:100), 2, function(v) {
  sample(c(rep(TRUE, 1750-nte), rep(FALSE, nte)), 1750, FALSE)  
}))
times <- c()
results <- matrix(NA, mc.reps, 1)
for (split.ind in 1:mc.reps) {
  t1 <- proc.time()
  trainYh2 <- trainYh + mvrnorm(n = 1500, mu = rep(0, 1331), Sigma = Sigma0)
  ttY <- rbind(trainYh2, testY)
  ##ttY <- rbind(trainY, testY)
  tr.inds <- train_splits[split.ind, ]
  tpred.new <- matrix(0, nte, nvox)
  fullbeta = matrix(0,nc=nvox,nr=ncol(c_ttF))
  for (vox.no in 1:nvox) {
    vox.ind <- usevox[vox.no]
    rule <- vox_prediction_rules[[vox.ind]]
    feat.inds <- which(rule$beta_coef != 0)
    if (length(feat.inds) > 0) {
      c_ttSub <- c_ttF[, feat.inds, drop = FALSE]
      teY <- ttY[!tr.inds, vox.ind]
      trY <- ttY[tr.inds, vox.ind]
      teX <- c_ttSub[!tr.inds, , drop = FALSE]
      trX <- c_ttSub[tr.inds, , drop = FALSE]
      bt <- ginv(trX) %*% trY
      Yh <- teX %*% bt
      fullbeta[feat.inds, vox.no] <- bt
      tpred.new[, vox.no] <- Yh
    }
  }
  voxcov = (t(ttY[tr.inds, usevox])%*%ttY[tr.inds, usevox])/(1750 - nte)
  voxfeatcov = (t(ttY[tr.inds, usevox]) %*% c_ttF[tr.inds, ])/(1750 - nte)
  voxCondCov = voxcov - voxfeatcov%*%fullbeta
  newCondVarS = (newCondVar + t(newCondVar))/2
  ridgeinv = solve(newCondVarS + diag(length(usevox))*3)
  score1 <- bestMatch(tpred.new, ttY[!tr.inds, ],'MSE',invCov = ridgeinv)
  results[split.ind, ] <- c(score1)
  (times[split.ind] <- (proc.time() - t1)["elapsed"])
}
results

