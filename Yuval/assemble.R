####
##  Reassemble files from pieces so we can save them on github
####

setwd("Yuval")
## v1_data.RData
p1 <- readRDS("v1_data_Rdata_1.rds")
p2 <- readRDS("v1_data_Rdata_2.rds")
p3 <- readRDS("v1_data_Rdata_3.rds")
p4 <- readRDS("v1_data_Rdata_4.rds")
p5 <- readRDS("v1_data_Rdata_5.rds")
fit_feat <- rbind(p1, p2, p3, p4, p5)
load("v1_data_Rdata_0.RData")
save(fit_feat, val_feat, SNRv1_corr, v1, file = "v1_data.RData")

## trainData.RData
load("v1_data.RData")
standards = apply(fit_feat,2,sd)
constF = which(standards==0)     
load('fit_trainS.RData') # loading fit_trainS
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
save(file = 'trainData.RData', c_trainF, trainY, c_testF, c_validF, testY, validY )

## bootstrap results
(fs <- paste0('temp_booting', paste(2:5), '.RData'))
load('temp_booting.RData')
rules_all <- rules
testpreds_all <-testpreds
rm(rules, testpreds)
for (f in fs) {
  load(f)
  if ('rules' %in% ls()) {
    rules_all <- c(rules_all, rules)
    rm(rules)
  }
  for (i in 1:1750) testpreds_all[[i]] <- rbind(testpreds_all[[i]],
                                                testpreds[[i]])
  rm(testpreds)
}
sapply(testpreds_all, nrow)
save('rules_all', 'testpreds_all', file = 'aggregate_booting.RData')

setwd("..")
