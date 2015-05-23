####
## Try EB-ident in V1
####

library(magrittr)
library(pracma)
library(MASS)
library(glmnet)
library(class)
source('transfer/source.R')
source('eb_ident/source.R')
source('eb_ident/eigenprism.R')
ddir <- "~/stat312data"
load(paste0(ddir, "/valid_index.RData"))
load(paste0(ddir, "/feature_valid.RData"))
dim(feature_valid) # 120 10921
load(paste0(ddir, "/valid_v1.RData"))
dim(valid_v1) # 1331 1560
temp <- read.csv(paste0(ddir, "/valid_resp_all.csv"), header = FALSE,
                 stringsAsFactors = FALSE)
temp %>% dim
temp %>% sapply(class)
filt <- temp %>% rowSums %>% is.na %>% `!`

Y <- temp[filt, ] %>% t %>% scale
dim(Y) # 120 22774
X <- feature_valid %>% scale
filt <- X %>% colSums %>% is.na %>% `!`
X <- X[, filt]
dim(X) # 120 10409

## Eigenprism procedure

res_EP <- eigenprisms(X, Y)
plot(res_EP$T2)
hist(res_EP$T2, breaks = 80)
