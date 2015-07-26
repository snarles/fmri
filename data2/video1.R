####
## Load 3500 training data
####

library(magrittr)
library(pracma, warn.conflicts = FALSE)
library(MASS)
library(glmnet, warn.conflicts = FALSE)
library(prodlim)
source('utils/zattach.R')
source('transfer/source.R')
source('eb_ident/eigenprism.R')
source('eb_ident/bayes_reg.R')
source('eb_ident/source.R')

gimage <- function(a) image(fliplr(t(a)), col = gray(0:20/20))

ddir <- '/home/snarles/stat312data/'
ddir <- '/home/rstudio/stat312data/'

list.files(ddir)
load(paste0(ddir, 'vidData1sec.RData'))
dim(vidTrainF) # 7128 8556
dim(vidTrainY) # 7128 2088

table(is.na(colSums(vidTrainY)))
plot(vidTrainY[, 1], type = "l")

cm <- cov(t(vidTrainY[1:2000, 20 * 1:100]))
dim(cm)

gimage(cm)

plot(vidTrainF[, 1], type = "l")
matplot(vidTrainF[, 100 * 1:5], type = "l")

fmeans <- colMeans(vidTrainF)

plot(fmeans, type = "l")




