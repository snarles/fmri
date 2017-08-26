source('Yuval/ident_setup.R')

tpred.new <- matrix(0, 1500, nvox)

nvox = 1250# up to 1250
usevox = order(SNRv1_corr,decreasing = TRUE)[1:nvox]

for (vox.no in 1:nvox) {
  vox.ind <- usevox[vox.no]
  rule <- vox_prediction_rules[[vox.ind]]
  feat.inds <- which(rule$beta_coef != 0)
  if (length(feat.inds) > 0) {
    trX <- c_trainF[, feat.inds, drop = FALSE]
    bt <- ginv(trX) %*% trainY[, vox.ind]
    Yh <- trX %*% bt
    tpred.new[, vox.no] <- Yh
  }
}

trainYh <- trainY; trainYh[, usevox] <- tpred.new
resid <- trainY - trainYh
sds <- apply(resid[, 1:1250], 2, sd)
(sd_m <- median(sds))
hist(sds)
