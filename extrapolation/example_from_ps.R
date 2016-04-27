## a bad test case

source("extrapolation/example2d.R")
source("extrapolation/mle_theory.R")
source("extrapolation/constrained_mle.R")

true_p_dist <- readRDS("extrapolation/testcase_ps.rds")

mean(true_p_dist^30)
mean(true_p_dist^90)

n <- 1000
Ys <- rbinom(n, 30, prob = sample(true_p_dist, n, TRUE))
mle_est <-  res_mixtools(Ys, 30)
pseq <- seq(0, 1, 1/300)
gu0 <- binned_gu(true_p_dist, pseq)
cm <- cons_mle_est(Ys, k, pseq, 0.01)
cmsamp <- sample(cm$ps, length(true_p_dist), replace = TRUE, prob = cm$gu)
plot(sort(true_p_dist), type = "l")
lines(sort(cmsamp), col = "red")

plot(gu0, type = "l")
lines(cm$gu, col = "red")
cm$of_gu(cm$gu)
cm$of_gu(gu0)


mean(true_p_dist^20)
mean(binmom(Ys, 30, 20))
est_moment(res, 20)
cm_est_moment(cm, 20)

mean(true_p_dist^30)
mean(binmom(Ys, 30, 30))
est_moment(res, 30)
cm_est_moment(cm, 30)

mean(true_p_dist^60)
est_moment(res, 60)
cm_est_moment(cm, 60)
