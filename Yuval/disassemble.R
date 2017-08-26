####
##  Cut files into pieces so we can save them on github
####

setwd("Yuval")

## v1_data.RData
load('v1_data.RData')
dim(fit_feat)
saveRDS(fit_feat[1:350, ], file = "v1_data_Rdata_1.rds")
saveRDS(fit_feat[351:700, ], file = "v1_data_Rdata_2.rds")
saveRDS(fit_feat[701:1050, ], file = "v1_data_Rdata_3.rds")
saveRDS(fit_feat[1051:1400, ], file = "v1_data_Rdata_4.rds")
saveRDS(fit_feat[1401:1750, ], file = "v1_data_Rdata_5.rds")
save(val_feat, SNRv1_corr, v1, file = "v1_data_Rdata_0.RData")

## trainData.RData
load('trainData.RData')
dim(c_trainF)
saveRDS(c_trainF[1:500, ], file = "trainData_RData_1.rds")
saveRDS(c_trainF[501:1000, ], file = "trainData_RData_2.rds")
saveRDS(c_trainF[1001:1500, ], file = "trainData_RData_3.rds")
saveRDS(trainY, file = "trainData_RData_4.rds")


