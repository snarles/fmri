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
p1 <- readRDS("trainData_RData_1.rds")
p2 <- readRDS("trainData_RData_2.rds")
p3 <- readRDS("trainData_RData_3.rds")
trainY <- readRDS("trainData_RData_4.rds")
c_trainF <- rbind(p1, p2, p3)
save(c_trainF, trainY, file = "trainData.RData")

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