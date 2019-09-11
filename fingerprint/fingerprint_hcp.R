library(pracma)
library(MASS)
library(class)
library(parallel)
library(lineId) ## use devtools::install('lineId')
source('extrapolation/kay_method2.R')
source('fingerprint/mi_est_pipeline.R')

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
nsubsample <- 20
all_inds <- matrix(nrow=nreps, ncol = nsubsample)

variable <- 'empirical_acc'
#variable <- 'mi_est'


for (i in 1:nreps) {
  # scores matrix
  set.seed(i)
  inds <- sample(nsubs, nsubsample)
  all_inds[i, ] <- inds
  # est_gsr_kl[i] <- mi_est_pipeline(gsr_kl[inds, inds])[[variable]]
  # est_ngsr_kl[i] <- mi_est_pipeline(ngsr_kl[inds, inds])[[variable]]
  # est_gsr_cor[i] <- mi_est_pipeline(gsr_cor[inds, inds])[[variable]]
  # est_ngsr_cor[i] <- mi_est_pipeline(ngsr_cor[inds, inds])[[variable]]
  est_gsr_cor[i] <- mi_est_pipeline2(gsr_cor[inds, inds])$mi_est
  est_ngsr_cor[i] <- mi_est_pipeline2(ngsr_cor[inds, inds])$mi_est
  acc_gsr_cor[i] <- mi_est_pipeline2(gsr_cor[inds, inds])$empirical_acc
  acc_ngsr_cor[i] <- mi_est_pipeline2(ngsr_cor[inds, inds])$empirical_acc
}

#results <- cbind(est_gsr_kl, est_ngsr_kl, est_gsr_cor, est_ngsr_cor)
#results <- cbind(acc_gsr_cor, acc_ngsr_cor, est_gsr_cor, est_ngsr_cor)
plot(jitter(acc_gsr_cor, 0.4), jitter(acc_ngsr_cor, 0.4), main = "Accuracy: GSR vs. no GSR")
abline(0, 1)

plot(est_gsr_cor, est_ngsr_cor, main = "I(X;Y): GSR vs. no GSR")
abline(0, 1)


plot(acc_gsr_cor - acc_ngsr_cor, est_gsr_cor - est_ngsr_cor, main = "delta I(X;Y) vs delta Acc")
abline(0, 0)
abline(v=0)



table(acc_gsr_cor - acc_ngsr_cor)
table(sign(est_gsr_cor - est_ngsr_cor))
table(acc_gsr_cor - acc_ngsr_cor, sign(est_gsr_cor - est_ngsr_cor ))
