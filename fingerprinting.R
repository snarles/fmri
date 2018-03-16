library(R.matlab)
library(pracma)
library(parallel)
mcc <- 20

fcs1 = readMat("~/Desktop/Results/All_Sub_REST1_TP.mat")
fcs2 = readMat("~/Desktop/Results/All_Sub_REST2_TP.mat")
fc = fcs1[[1]][1,,]
hist(fc[upper.tri(fc)])
upper.tri(matrix(1:4, 2, 2))

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

t1 <- proc.time()
nsb <- 20
klds <- matrix(0, nsb, nsb)
for (i in 1:nsb) {
  for (j in 1:nsb) {
    vstuff <- solve(mats2[i,,], mats1[j,,])
    klds[i,j] <- sum(diag(vstuff)) - log(det(vstuff))
  }
}

mean(apply(-klds, 1, which.max) == 1:nrow(klds))
proc.time() - t1

ijs <- apply(cbind(rep(1:nsubs, each=nsubs), rep(1:nsubs, nsubs)), 1, list)
get_kl_subroutine <- function(ij) {
  i <- ij[[1]][1]
  j <- ij[[1]][2]
  vstuff <- solve(mats1[i,,], mats2[j,,])
  sum(diag(vstuff)) - log(det(vstuff))
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
df_obs <- 900
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


nboot <- 20
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
bt <- nnls::nnls(Xmat, accs_sub)
accs_cvr <- Xpred %*% bt$x

library(minpack.lm)

regr <- nlsLM(accs ~ c + b * exp(-tt/x), data = list(accs = accs_sub, tt=1:nsb), start = list(c = 0.5, b = 0.1, x = 70))
regr
accs_E <- coef(regr)['c'] + coef(regr)['b'] * exp(-(1:nsubs2)/coef(regr)['x'])

layout(1)
plot(accs, type = "l", ylim = c(0,1), lwd =2)
lines(accs_par2, col = "green", lwd =2)
lines(accs_YB, col = "blue", lwd =2)
lines(accs_W, col = "purple", lwd =2)
lines(accs_KK, col = "orange", lwd =2)
lines(accs_E, col = "pink", lwd =2)
lines(accs_cvr, col = "brown", lwd = 2)
lines(accs_sub, col = "red", lwd = 4)

#saveRDS(r12, 'par2example.rds')
