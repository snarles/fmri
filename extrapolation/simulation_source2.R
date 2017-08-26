## Compare original method for Ku extrapolation and top-m rank method (zero-one loss)

source('extrapolation/simulation_models.R')
library(lineId)

coefficient_matrix_type1 <- function(k1, k2 = k1, d) {
  ans <- matrix(0, k2 - k1 + 1, d + 1)
  ii <- row(ans); jj <- col(ans)
  ans <- matrix((k1 + ii - 2)/(k1 + ii + jj - 3), nrow = nrow(ans))
  ans[ii == 1 & jj == 1] <- 1 
  ans
}

coefficient_matrix_type2 <- function(k, d) {
  ans <- matrix(0, k, d + 1)
  ii <- row(ans); jj <- col(ans)
  ans <- matrix(beta(k-1-(ii-1)+(jj-1), ii)/beta(k-1-(ii-1), ii), nrow = nrow(ans))
  ans[ii == nrow(ans) & jj == 1] <- 1
  ans
}

qda_params1 <- list(
  k = 200,
  r_train = 40,
  r_test = 1,
  p = 2,
  wishart_par = 3.5,
  snr = 30
)
qda_params2 <- modifyList(qda_params1, list(k = 1000))

## Timing data
## k= 1000 takes 1.44 seconds, then scales as ^2
## k= 4000 takes 14.8s
# t1 <- proc.time()
# res <- do.call(simulate_qda, qda_params2)
# mc_rates <- resample_misclassification(res$plikes, res$i_chosen, 1:ncol(res$plikes))
# proc.time() - t1
# rm(res)

mc_ratess <- sapply(1:100, 
                    function(ind) {
                      res <- do.call(simulate_qda, qda_params2)
                      resample_misclassification(res$plikes, res$i_chosen, 1:ncol(res$plikes))
                    })
dim(mc_ratess)
matplot(mc_ratess, type = "l", ylim = c(0, 1))
mc_rates <- rowMeans(mc_ratess)
sds <- apply(mc_ratess, 1, sd)
plot(sds, type = "l")

## extrapolate from fewer classes

mc.reps <- 100
errs_tot_s <- numeric()
errs_top_s <- numeric()

for (ind in 1:mc.reps) {
  out <- data.frame(true = mc_rates)
  
  for (d in 1:20) {
    pb <- coefficient_matrix_type1(1, qda_params2$k, d = d)
    dmults <- c(4, 5, 6)
    res <- do.call(simulate_qda, qda_params1)
    errs1 <- resample_misclassification(res$plikes, res$i_chosen, 1:qda_params1$k)
    errs2 <- 1 - topk_scores(res$plikes, res$i_chosen, qda_params1$k)/qda_params1$k
    mat1 <- coefficient_matrix_type1(1, qda_params1$k, d)
    mat2 <- coefficient_matrix_type2(qda_params1$k, d)
    # mat1 <- mat1[nrow(mat1):1, ]
    # errs1 <- rev(errs1)
    # mat1 <- mat1[1:190, ]
    # errs1 <- errs1[1:190]
    res1 <- as.numeric(pb %*% pinv(mat1) %*% errs1)
    nm1 <- paste0("m1_", d)
    out[[nm1]] <- res1
    
    for (dm in dmults) {
      m2 <- mat2[1:(d * dm), ]
      e2 <- errs2[1:(d * dm)]
      res2 <- as.numeric(pb %*% pinv(m2) %*% e2)
      nm <- paste0("m2_", d, "_dm_", dm)
      out[[nm]] <- res2
    }
  }
  
  errs_top <- sapply(out, function(v) sum(abs(v - mc_rates)[-(1:qda_params1$k)]), 
                     USE.NAMES = TRUE)
  sort(errs_top)
  errs_top_s <- rbind(errs_top_s, errs_top)
  errs_tot <- sapply(out, function(v) sum(abs(v - mc_rates)), USE.NAMES = TRUE)
  sort(errs_tot)
  errs_tot_s <- rbind(errs_tot_s, errs_tot)
}
errs_tot_av <- colMeans(errs_tot_s)
sort(errs_tot_av)
errs_top_av <- colMeans(errs_top_s)
sort(errs_top_av)

# matplot(out[names(sort(errs_tot_av))[1:5]], ylim = c(0, 1), type = "l")
# legend(400, 0.6, col = 1:5, lty = 1:5, legend = names(sort(errs_tot_av)[1:5]))

