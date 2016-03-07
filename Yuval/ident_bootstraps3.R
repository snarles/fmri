####
##  Charles' stuff using Yuval's stuff
##  Bootstrap the accuracy numbers using the same training data
####

source('Yuval/ident_setup.R')

###
## PARAMS
###

nvox = 100# up to 1250
nsubsample <- 1250
nboot <- 20
ms <- c(50, 100, 150, 200, 250)


###
## Setup:: estimate inverse cov mat
###
nm <- length(ms)
ihats <- matrix(NA, nboot, nm)

usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]
testpred = getPreds(c_testF,ind_struct=vox_prediction_rules,voxind=usevox)
newCondVar =genCondVar(1:10409,usevox)
newCondVarS = (newCondVar + t(newCondVar))/2
ridgeinv = solve(newCondVarS + diag(length(usevox))*3)

###
##  For loop involving nboot
###
for (boot.ind in 1:nboot) {
  inds_sub <- sample(1500, nsubsample, FALSE)
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
  ##(score1 <- bestMatch(tpred.new, testY, 'MSE',invCov = ridgeinv))
  scores <- matchResults(tpred.new, testY, 'MSE',invCov = ridgeinv, rankit = FALSE)
  temp <- t(apply(scores, 1, rank))
  sum(diag(temp) == 1)/250
  
  for (ind in 1:nm) {
    m <- ms[ind]
    (p <- resample_misclassification(-scores, 1:250, m, replace = FALSE))
    (p_prime <- rbinom(1, 250, p)/250)
    (IH <- Ihat_LI(p_prime, m1, 20))
    ihats[boot.ind, ind] <- IH
  }
}



###
## SAVE RESULTS
###

fname <- paste0('Yuval/ident_bootstraps3_nvox', nvox, '.rds')
saveRDS(ihats, file = fname)

###
## ANALYZE RESULTS
###
# 
# m1 <- 200
# m2 <- 50
# (p <- resample_misclassification(-scores, 1:250, m1, replace = FALSE))
# (p_prime <- rbinom(1, 250, p)/250)
# (IH <- Ihat_LI(p_prime, m1, 10))
# piK(sqrt(2 * IH), m2)
# (resample_misclassification(-scores, 1:250, m2, replace = FALSE))
# 
