library(R.matlab)
library(pracma)
library(parallel)
mcc <- 20

#fcs1 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST1_RP.mat")
#fcs2 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST2_RP.mat")

#fcs1 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST1_errts.mat")
#fcs2 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST2_errts.mat")

fcs1 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_REST1_errts.mat")
fcs2 = readMat("/data/HCP_preproc/fingerprinting/Results/All_Sub_SOCIAL_errts.mat")
dim(fcs2[[1]])

#fcs1 = readMat("~/Desktop/Results/All_Sub_REST1_TP.mat")
#fcs2 = readMat("~/Desktop/Results/All_Sub_REST2_TP.mat")

#fcs1 = readMat("~/Desktop/Results/All_Sub_Rest1.mat")
#fcs2 = readMat("~/Desktop/Results/All_Sub_Rest2.mat")

#fc = fcs1[[1]][1,,]
#hist(fc[upper.tri(fc)])
#upper.tri(matrix(1:4, 2, 2))

fc_flatten <- function(arr) {
  apply(arr, 1, function(a) {
    a[upper.tri(a)]
  })
}

nsubs <- nrow(fcs1[[1]])
## shuffle subjects
shufind <- sample(nsubs, nsubs, replace = FALSE)

fcz1 <- fc_flatten(fcs1[[1]])[,shufind]
fcz2 <- fc_flatten(fcs2[[1]])[,shufind]

#r12 <- atanh(cor(fcz1, fcz2))

r12 <- cor(fcz1, fcz2)
#r12 <- -pdist2(t(fcz1), t(fcz2))/100 + 1
image(r12)


display_results <- function(r12) {
  print(mean(apply(r12, 1, which.max) == 1:nrow(r12)))
  layout(t(t(1:2)))
  plot(density(diag(r12), bw = 'SJ'), xlim = c(min(r12), max(r12)), ylim = c(0, 12/(max(r12) - min(r12))), lwd = 2)
  for (i in 1:33 * 10) {
    lines(density(r12[i, -i], bw = 'SJ'), col = "red")
  }
  lines(density(r12[upper.tri(r12)]), col = "blue", lwd = 2)
  
  r_mean <- rowSums(r12 - diag(diag(r12)))/(ncol(r12) - 1)
  r_std <- apply((r12 - diag(diag(r12))), 2, std)
  r_sub2 <- (r12 - r_mean)/r_std
  plot(density(diag(r_sub2), bw = 'SJ'), xlim = c(min(r_sub2), max(r_sub2)), ylim = c(0, 12/(max(r_sub2) - min(r_sub2))), lwd = 2)
  for (i in 1:33 * 10) {
    lines(density(r_sub2[i, -i], bw = 'SJ'), col = "red")
  }
  lines(density(r_sub2[upper.tri(r12)], bw = 'SJ'), col = "blue", lwd = 2)

}

display_results(r12)

####
##  Extrapolation
####

source("par2/objective_function.R")

library(lineId)

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/kay_method.R")
nsb <- 100

accs <- 1 - resample_misclassification(r12, 1:nsubs, 1:nsubs)
accs_sub <- 1 - resample_misclassification(r12[1:nsb, 1:nsb], 1:nsb, 1:nsb)
#accs_sub <- accs[1:nsb]
accs_par2 <- par2_extrapolate(1:nsb, accs_sub, 1:nsubs)
r_mean <- rowSums(r12 - diag(diag(r12)))/(ncol(r12) - 1)
r_std <- apply((r12 - diag(diag(r12))), 2, std)
r_sub2 <- (r12 - r_mean)/r_std
muh_YB <- mean(diag(r_sub2))
tau_YB <- var(diag(r_sub2))
accs_YB <- par2_acc_k(1:nsubs, muh_YB, tau_YB)

mean_id <- mean(diag(r12))
mean_nonid <- mean(r12[upper.tri(r12)])
tau_W <- var(c(diag(r12) - mean_id, r12[upper.tri(r12)] - mean_nonid))
accs_W <- par2_acc_k(1:nsubs, (mean_id - mean_nonid)/sqrt(tau_W), 1)

accs_KK <- kernel_extrap(r12[1:nsb, 1:nsb], 1:nsubs, bw = 'ucv')

library(minpack.lm)

regr <- nlsLM(accs ~ c + b * exp(-tt/x), data = list(accs = accs_sub, tt=1:nsb), start = list(c = 0.5, b = 0.1, x = 70))
regr
accs_E <- coef(regr)['c'] + coef(regr)['b'] * exp(-(1:nsubs)/coef(regr)['x'])

layout(1)
plot(accs, type = "l", ylim = c(0,1), lwd =2)
lines(accs_par2, col = "green", lwd =2)
lines(accs_YB, col = "blue", lwd =2)
lines(accs_W, col = "purple", lwd =2)
lines(accs_KK, col = "orange", lwd =2)
lines(accs_E, col = "pink", lwd =2)
lines(accs_sub, col = "red", lwd = 4)


####
## Try Kld
####

mats1 <- fcs1[[1]]
mats2 <- fcs2[[1]]

ijs <- apply(cbind(rep(1:nsubs, each=nsubs), rep(1:nsubs, nsubs)), 1, list)
get_kl_subroutine <- function(ij) {
  i <- ij[[1]][1]
  j <- ij[[1]][2]
  vstuff <- solve(mats1[i,,], mats2[j,,])
  vstuff2 <- solve(mats2[i,,], mats1[j,,])  
  pmin(sum(diag(vstuff)) - log(det(vstuff)), sum(diag(vstuff2)) - log(det(vstuff2)))
}

t1 <- proc.time()
kls <- mclapply(ijs, get_kl_subroutine, mc.cores = mcc)
proc.time() - t1

klds <- t(matrix(unlist(kls), nsubs, nsubs))
display_results(-klds)

####
## Simulate data
####

library(MCMCpack)
nrois <- 268
nsubs2 <- 1000
ncomps <- 20
df_sig <- 300
df_obs <- 2100 # acc 0.855
#df_obs <- 1800 # acc 0.795
#df_obs <- 1500 # acc 0.735
#df_obs <- 1200 # acc .57
#df_obs <- 900 # acc .5
#df_obs <- 600 # acc .3
diralpha <- rep(.2, ncomps)

templates <- array(0, dim = c(ncomps, nrois, nrois))
for (i in 1:ncomps) {
  templates[i,,] <- rWishart(1, df_sig, eye(nrois))/df_sig
}
signals <- array(0, dim = c(nsubs2, nrois, nrois))

wts <- rdirichlet(nsubs2, diralpha)
dim(wts)

for (j in 1:ncomps) {
  for (i in 1:nsubs2) {
    signals[i,,] <- signals[i,,] + wts[i,j] * templates[j,,]
  }
}

obs1 <- array(t(apply(signals, 1, function(v) cov2cor(rWishart(1, df_obs, v)[,,1]/df_obs))), dim = c(nsubs2, nrois, nrois))
obs2 <- array(t(apply(signals, 1, function(v) cov2cor(rWishart(1, df_obs, v)[,,1]/df_obs))), dim = c(nsubs2, nrois, nrois))

#plot(obs1[1,,], signals[1,,], pch = ".", xlim = c(-.2, .2), ylim = c(-.2,.2))
#plot(obs1[1,,], signals[2,,], pch = ".", xlim = c(-.2, .2), ylim = c(-.2,.2))

obs1f <- fc_flatten(obs1)
obs2f <- fc_flatten(obs2)
r12 <- cor(obs1f, obs2f)

display_results(r12)

####
## Extrapolation setup
####

source("par2/objective_function.R")

library(lineId)

library(minpack.lm)

source("approximation/gaussian_identity_finsam.R")
source("approximation/gaussian_identity_finsam2.R")
source("extrapolation/kay_method.R")
source("extrapolation/ku_source.R")
source("extrapolation/basis_source.R")


nsb <- 200
kref <- 1:nsb
Ktarg <- 1:nsubs2
(max.mu <- (qnorm(1- 1/(max(kref)^2))))
lsub2 <- length(kref)/2
(bdwids <- seq(0.1, 1, 0.1))
basis_sets <- lapply(bdwids, function(bd) get_basis_mat(max.mu, bd))

sub_basis_sets <- lapply(basis_sets, function(set1) {
  list(Xmat = set1$Xmat[1:lsub2, ], Xtarg = set1$Xmat[length(kref), , drop = FALSE])
})
nboot <- 20

####
##  Extrapolation
####

accs <- 1 - resample_misclassification(r12, 1:nsubs2, 1:nsubs2)
accs_sub <- 1 - resample_misclassification(r12[1:nsb, 1:nsb], 1:nsb, 1:nsb)
#accs_sub <- accs[1:nsb]
accs_par2 <- par2_extrapolate(1:nsb, accs_sub, 1:nsubs2, verbose = TRUE)
r_mean <- rowSums(r12 - diag(diag(r12)))/(ncol(r12) - 1)
r_std <- apply((r12 - diag(diag(r12))), 2, std)
r_sub2 <- (r12 - r_mean)/r_std
muh_YB <- mean(diag(r_sub2))
tau_YB <- var(diag(r_sub2))
accs_YB <- par2_acc_k(1:nsubs2, muh_YB, tau_YB)

mean_id <- mean(diag(r12))
mean_nonid <- mean(r12[upper.tri(r12)])
tau_W <- var(c(diag(r12) - mean_id, r12[upper.tri(r12)] - mean_nonid))
accs_W <- par2_acc_k(1:nsubs2, (mean_id - mean_nonid)/sqrt(tau_W), 1)

accs_KK <- kernel_extrap(r12[1:nsb, 1:nsb], 1:nsubs2, bw = 'ucv')



boot_accs <- matrix(NA, nboot, lsub2)
for (ii in 1:nboot) {
  subinds <- sample(nsb, lsub2, replace = FALSE)
  boot_accs[ii, ] <- 1- resample_misclassification(r12[subinds, subinds], 1:lsub2, 1:lsub2)
}

all_sub_preds <- t(apply(boot_accs, 1, bdwid_all_preds, basis_sets = sub_basis_sets))
cv_curve <- sqrt(colMeans((all_sub_preds - accs_sub[length(kref)])^2))
plot(bdwids, cv_curve, type ="l")
sel_ind <- which.min(cv_curve)
Xmat <- basis_sets[[sel_ind]]$Xmat
Xpred <- basis_sets[[sel_ind]]$Xtarg
pen <- 1000
bt <- nnls::nnls(rbind(Xmat, rep(pen, ncol(Xmat))), c(accs_sub, pen))
sum(bt$x)
accs_cvr <- Xpred %*% bt$x

regr <- nlsLM(accs ~ c + b * exp(-tt/x), data = list(accs = accs_sub, tt=1:nsb), start = list(c = 0.5, b = 0.1, x = 70))
accs_E <- coef(regr)['c'] + coef(regr)['b'] * exp(-(1:nsubs2)/coef(regr)['x'])

#regr <- nlsLM(accs ~ c0 + c1 * tt + b * exp(-tt/x), data = list(accs = accs_sub, tt=1:nsb), start = list(c0 = 0.5, c1 = 0, b = 0.1, x = 70))
#accs_E2 <- coef(regr)['c0'] + coef(regr)['c1'] * (1:nsubs2) + coef(regr)['b'] * exp(-(1:nsubs2)/coef(regr)['x'])

layout(1)
plot(accs, type = "l", ylim = c(0,1), lwd =2)
lines(accs_par2, col = "green", lwd =2)
lines(accs_YB, col = "blue", lwd =2)
lines(accs_W, col = "purple", lwd =2)
lines(accs_KK, col = "orange", lwd =2)
lines(accs_E, col = "pink", lwd =2)
#lines(accs_E2, col = "blue", lwd =2)
lines(accs_cvr, col = "brown", lwd = 2)
lines(accs_sub, col = "red", lwd = 4)

#saveRDS(r12, 'par2example.rds')


results <- rbind(true_accs = accs, sub_accs = c(accs_sub, rep(NA, nsubs2-nsb)), 
      par2 = accs_par2, yb = accs_YB, wl = accs_W, kay = accs_KK, exex = accs_E, 
      cvr = accs_cvr[, 1])


####
##  Run many experiments
####

simulate_and_run <- function(rseed, ncomps, diralpha, df_sig, df_obs, nboot = 20) {
  set.seed(rseed)
  templates <- array(0, dim = c(ncomps, nrois, nrois))
  for (i in 1:ncomps) {
    templates[i,,] <- rWishart(1, df_sig, eye(nrois))/df_sig
  }
  signals <- array(0, dim = c(nsubs2, nrois, nrois))
  
  wts <- rdirichlet(nsubs2, diralpha)
  dim(wts)
  
  for (j in 1:ncomps) {
    for (i in 1:nsubs2) {
      signals[i,,] <- signals[i,,] + wts[i,j] * templates[j,,]
    }
  }
  
  obs1 <- array(t(apply(signals, 1, function(v) cov2cor(rWishart(1, df_obs, v)[,,1]/df_obs))), dim = c(nsubs2, nrois, nrois))
  obs2 <- array(t(apply(signals, 1, function(v) cov2cor(rWishart(1, df_obs, v)[,,1]/df_obs))), dim = c(nsubs2, nrois, nrois))
  
  obs1f <- fc_flatten(obs1)
  obs2f <- fc_flatten(obs2)
  r12 <- cor(obs1f, obs2f)
  
  
  ## extrapolate
  accs <- 1 - resample_misclassification(r12, 1:nsubs2, 1:nsubs2)
  accs_sub <- 1 - resample_misclassification(r12[1:nsb, 1:nsb], 1:nsb, 1:nsb)
  #accs_sub <- accs[1:nsb]
  accs_par2 <- par2_extrapolate(1:nsb, accs_sub, 1:nsubs2, verbose = TRUE)
  r_mean <- rowSums(r12 - diag(diag(r12)))/(ncol(r12) - 1)
  r_std <- apply((r12 - diag(diag(r12))), 2, std)
  r_sub2 <- (r12 - r_mean)/r_std
  muh_YB <- mean(diag(r_sub2))
  tau_YB <- var(diag(r_sub2))
  accs_YB <- par2_acc_k(1:nsubs2, muh_YB, tau_YB)
  
  mean_id <- mean(diag(r12))
  mean_nonid <- mean(r12[upper.tri(r12)])
  tau_W <- var(c(diag(r12) - mean_id, r12[upper.tri(r12)] - mean_nonid))
  accs_W <- par2_acc_k(1:nsubs2, (mean_id - mean_nonid)/sqrt(tau_W), 1)
  
  accs_KK <- kernel_extrap(r12[1:nsb, 1:nsb], 1:nsubs2, bw = 'ucv')
  
  
  
  boot_accs <- matrix(NA, nboot, lsub2)
  for (ii in 1:nboot) {
    subinds <- sample(nsb, lsub2, replace = FALSE)
    boot_accs[ii, ] <- 1- resample_misclassification(r12[subinds, subinds], 1:lsub2, 1:lsub2)
  }
  
  all_sub_preds <- t(apply(boot_accs, 1, bdwid_all_preds, basis_sets = sub_basis_sets))
  cv_curve <- sqrt(colMeans((all_sub_preds - accs_sub[length(kref)])^2))
  #plot(bdwids, cv_curve, type ="l")
  sel_ind <- which.min(cv_curve)
  Xmat <- basis_sets[[sel_ind]]$Xmat
  Xpred <- basis_sets[[sel_ind]]$Xtarg
  pen <- 1000
  bt <- nnls::nnls(rbind(Xmat, rep(pen, ncol(Xmat))), c(accs_sub, pen))
  sum(bt$x)
  accs_cvr <- Xpred %*% bt$x
  regr <- nlsLM(accs ~ c + b * exp(-tt/x), data = list(accs = accs_sub, tt=1:nsb), start = list(c = 0.5, b = 0.1, x = 70))
  accs_E <- coef(regr)['c'] + coef(regr)['b'] * exp(-(1:nsubs2)/coef(regr)['x'])
  results <- rbind(true_accs = accs, sub_accs = c(accs_sub, rep(NA, nsubs2-nsb)), 
                   par2 = accs_par2, yb = accs_YB, wl = accs_W, kay = accs_KK, exex = accs_E, 
                   cvr = accs_cvr[, 1])
  return(results)
}





results <- simulate_and_run(0, ncomps, diralpha, 300, df_obs)



library(parallel)
mcc <- 20
mcc

results300 <- mclapply(1:mcc, function(i) simulate_and_run(i, ncomps, diralpha, df_sig, 300), mc.cores = mcc)
results600 <- mclapply(1:mcc, function(i) simulate_and_run(i, ncomps, diralpha, df_sig, 600), mc.cores = mcc)
results900 <- mclapply(1:mcc, function(i) simulate_and_run(i, ncomps, diralpha, df_sig, 900), mc.cores = mcc)
results1200 <- mclapply(1:mcc, function(i) simulate_and_run(i, ncomps, diralpha, df_sig, 1200), mc.cores = mcc)
results1500 <- mclapply(1:mcc, function(i) simulate_and_run(i, ncomps, diralpha, df_sig, 1500), mc.cores = mcc)

# plotresults <- function(results) {
#   accs <- results["true_accs", ]
#   accs_sub <- results["sub_accs", 1:nsb]
#   accs_par2 <- results["par2", ]
#   accs_YB <- results["yb", ]
#   accs_W <- results["wl", ]
#   accs_KK <- results["kay", ]
#   accs_E <- results["exex", ]
#   accs_cvr <- results["cvr", ]
#   layout(1)
#   plot(accs, type = "l", ylim = c(0,1), lwd =2)
#   lines(accs_par2, col = "green", lwd =2)
#   lines(accs_YB, col = "blue", lwd =2)
#   lines(accs_W, col = "purple", lwd =2)
#   lines(accs_KK, col = "orange", lwd =2)
#   lines(accs_E, col = "pink", lwd =2)
#   #lines(accs_E2, col = "blue", lwd =2)
#   lines(accs_cvr, col = "brown", lwd = 2)
#   lines(accs_sub, col = "red", lwd = 4)
# }

plotresults <- function(results) {
  accs <- results["true_accs", ]
  accs_sub <- results["sub_accs", 1:nsb]
  accs_par2 <- results["par2", ]
  accs_YB <- results["yb", ]
  accs_W <- results["wl", ]
  accs_KK <- results["kay", ]
  accs_E <- results["exex", ]
  accs_cvr <- results["cvr", ]
  layout(1)
  plot(accs, type = "l", ylim = c(0,1), lwd =2)
  lines(accs_par2, col = "green", lwd =2)
  #lines(accs_YB, col = "blue", lwd =2)
  lines(accs_W, col = "purple", lwd =2)
  #lines(accs_KK, col = "orange", lwd =2)
  lines(accs_E, col = "orange", lwd =2)
  #lines(accs_E2, col = "blue", lwd =2)
  #lines(accs_cvr, col = "brown", lwd = 2)
  lines(accs_sub, col = "red", lwd = 4)
  legend(0, 0.4, legend = c("true", "subset", "waller", "expo", "par2"), lwd = 2, 
         col = c("black", "red", "purple", "orange", "green"))
}

summaryresults <- function(results, sel_rows = 1:nsubs2) {
  diffs <- t(t(results[3:8, ]) - results[1, ])
  sqrt(rowMeans(diffs[, sel_rows, drop = FALSE]^2))
}


i <- 0
i <- i + 1; plotresults(results1500[[i]]); title(paste("High SNR", i))

i <- 0
i <- i + 1; plotresults(results1200[[i]]); title(paste("Med. High SNR", i))

i <- 0
i <- i + 1; plotresults(results900[[i]]); title(paste("Medium SNR", i))

i <- 0
i <- i + 1; plotresults(results600[[i]]); title(paste("Low SNR", i))


inds <- c(3, 5, 1)
rowMeans(sapply(results300, summaryresults))[inds]
rowMeans(sapply(results600, summaryresults))[inds]
rowMeans(sapply(results900, summaryresults))[inds]
rowMeans(sapply(results1200, summaryresults))[inds]
rowMeans(sapply(results1500, summaryresults))[inds]

# > rowMeans(sapply(results300, summaryresults))
# par2         yb         wl        kay       exex        cvr 
# 0.04563530 0.53901913 0.51577173 0.17302000 0.09060635 0.04282269 
# > rowMeans(sapply(results600, summaryresults))
# par2         yb         wl        kay       exex        cvr 
# 0.05083113 0.42824756 0.37413371 0.14253290 0.05973530 0.05224761 
# > rowMeans(sapply(results900, summaryresults))
# par2         yb         wl        kay       exex        cvr 
# 0.04284083 0.32205314 0.25083351 0.10488390 0.05622053 0.06699989 
# > rowMeans(sapply(results1200, summaryresults))
# par2         yb         wl        kay       exex        cvr 
# 0.05124250 0.23971899 0.15795199 0.07163714 0.05269437 0.06016347 
# > rowMeans(sapply(results1500, summaryresults))
# par2         yb         wl        kay       exex        cvr 
# 0.04383392 0.17575085 0.08736941 0.04724770 0.04430104 0.05214578


rowMeans(sapply(results300, summaryresults, sel_rows = 1000))
rowMeans(sapply(results600, summaryresults, sel_rows = 1000))
rowMeans(sapply(results900, summaryresults, sel_rows = 1000))
rowMeans(sapply(results1200, summaryresults, sel_rows = 1000))
rowMeans(sapply(results1500, summaryresults, sel_rows = 1000))

# > rowMeans(sapply(results300, summaryresults, sel_rows = 1000))
# par2         yb         wl        kay       exex        cvr 
# 0.04645720 0.67675421 0.65123713 0.25369290 0.14856859 0.04105121 
# > rowMeans(sapply(results600, summaryresults, sel_rows = 1000))
# par2         yb         wl        kay       exex        cvr 
# 0.06710622 0.57053176 0.50502455 0.22507151 0.09167060 0.06306193 
# > rowMeans(sapply(results900, summaryresults, sel_rows = 1000))
# par2         yb         wl        kay       exex        cvr 
# 0.05479273 0.44797352 0.35848232 0.17523901 0.08775081 0.09510723 
# > rowMeans(sapply(results1200, summaryresults, sel_rows = 1000))
# par2         yb         wl        kay       exex        cvr 
# 0.07262386 0.34477589 0.23945394 0.12552560 0.07998243 0.09107651 
# > rowMeans(sapply(results1500, summaryresults, sel_rows = 1000))
# par2         yb         wl        kay       exex        cvr 
# 0.06899284 0.25979021 0.14339937 0.08337639 0.07387766 0.08173827 