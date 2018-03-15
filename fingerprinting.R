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
  print(mean(apply(r12, 1, which.max) == 1:339))
  layout(t(t(1:2)))
  plot(density(diag(r12)), xlim = c(min(r12), max(r12)), ylim = c(0, 12/(max(r12) - min(r12))), lwd = 2)
  for (i in 1:33 * 10) {
    lines(density(r12[i, -i]), col = "red")
  }
  lines(density(r12[upper.tri(r12)]), col = "blue", lwd = 2)
  
  r_mean <- rowSums(r12 - diag(diag(r12)))/(ncol(r12) - 1)
  r_std <- apply((r12 - diag(diag(r12))), 2, std)
  r_sub2 <- (r12 - r_mean)/r_std / 15
  plot(density(diag(r_sub2)), xlim = c(min(r_sub2), max(r_sub2)), ylim = c(0, 12/(max(r_sub2) - min(r_sub2))), lwd = 2)
  for (i in 1:33 * 10) {
    lines(density(r_sub2[i, -i], bw = 0.03), col = "red")
  }
  lines(density(r_sub2[upper.tri(r12)]), col = "blue", lwd = 2)

}

display_results(r12)

####
##  Extrapolation
####

source("par2/objective_function.R")

library(lineId)

help("resample_misclassification")

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
lines(accs_sub, col = "red", lwd = 4)
lines(accs_par2, col = "green", lwd =2)
lines(accs_YB, col = "blue", lwd =2)
lines(accs_W, col = "purple", lwd =2)
lines(accs_KK, col = "orange", lwd =2)
lines(accs_E, col = "pink", lwd =2)


####
## Simulate data
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
