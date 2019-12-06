library(pracma)
library(MASS)
library(class)
library(parallel)
library(lineId) ## use devtools::install('lineId')
library(R.matlab)
library(FNN)

source('extrapolation/kay_method2.R')
source('fingerprint/mi_est_pipeline.R')


fc_flatten <- function(arr) {
  apply(arr, 1, function(a) {
    a[upper.tri(a)]
  })
}


fcs1gsr = fc_flatten(readMat(gzfile("~/Documents/fingerprint_data/FL_Matrices/All_Sub_REST1_FL_GSR.mat.gz"))[[1]])
fcs2gsr = fc_flatten(readMat(gzfile("~/Documents/fingerprint_data/FL_Matrices/All_Sub_REST2_FL_GSR.mat.gz"))[[1]])
fcs1ngsr = fc_flatten(readMat(gzfile("~/Documents/fingerprint_data/FL_Matrices/All_Sub_REST1_FL.mat.gz"))[[1]])
fcs2ngsr = fc_flatten(readMat(gzfile("~/Documents/fingerprint_data/FL_Matrices/All_Sub_REST2_FL.mat.gz"))[[1]])

mis_gsr <- sapply(1:dim(fcs1gsr)[1], function(i) mutinfo(fcs1gsr[i, ], fcs2gsr[i, ]))
mis_ngsr <- sapply(1:dim(fcs1gsr)[1], function(i) mutinfo(fcs1ngsr[i, ], fcs2ngsr[i, ]))


gsr_kl <- -as.matrix(read.csv('fingerprint/FL_GSR_kl.csv'))
gsr_cor <- as.matrix(read.csv('fingerprint/FL_GSR_cor.csv'))
ngsr_kl <- -as.matrix(read.csv('fingerprint/FL_noGSR_kl.csv'))
ngsr_cor <- as.matrix(read.csv('fingerprint/FL_noGSR_cor.csv'))


nsubs <- dim(gsr_kl)[1] # 338






est_gsr_kl <- numeric()
est_ngsr_kl <- numeric()
est_gsr_cor <- numeric()
est_ngsr_cor <- numeric()
acc_gsr_cor <- numeric()
acc_ngsr_cor <- numeric()


nreps <- 100
nsubsample <- 40
all_inds <- matrix(nrow=nreps, ncol = nsubsample)

variable <- 'empirical_acc'
#variable <- 'mi_est'

mi_est_pipeline2(gsr_kl)$empirical_acc
mi_est_pipeline2(ngsr_kl)$empirical_acc
tot_mi_gsr <- mi_est_pipeline2(gsr_kl)$mi_est
tot_mi_ngsr <- mi_est_pipeline2(ngsr_kl)$mi_est


##
d_gsr <- tot_mi_gsr/sum(mis_gsr) * dim(fcs1gsr)[1]
d_ngsr <- tot_mi_ngsr/sum(mis_ngsr) * dim(fcs1gsr)[1]

convert_xx_to_xy <- function(mi_xx, d_est) {
  mi_per <- mi_xx/d_est
  rho_est <- (1-exp(-2 * mi_per))^0.25
  mi_per_xy <- -.5 * log(1-rho_est^2)
  d_est * mi_per_xy
}

mi_xy_gsr <- convert_xx_to_xy(tot_mi_gsr, d_gsr)
mi_xy_ngsr <- convert_xx_to_xy(tot_mi_ngsr, d_ngsr)


data.frame(tot_mi_gsr, tot_mi_ngsr, mi_xy_gsr, mi_xy_ngsr, d_gsr, d_ngsr)

# tot_mi_gsr tot_mi_ngsr mi_xy_gsr mi_xy_ngsr    d_gsr  d_ngsr
# 1   11.41317    12.02459  26.95941   24.18459 75.27071 50.7938