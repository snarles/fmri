## Compare original method for Ku extrapolation and top-m rank method (zero-one loss)

source('extrapolation/simulation_models.R')

coefficient_matrix_type1 <- function(k1, k2, d) {
  ans <- matrix(0, k2 - k1 + 1, d + 1)
  ii <- row(ans); jj <- col(ans)
  ans <- matrix((k1 + ii - 2)/(k1 + ii + jj - 3), nrow = nrow(ans))
  ans
}

coefficient_matrix_type2 <- function(k, d) {
  ans <- matrix(0, k - 1, d + 1)
  ii <- row(ans); jj <- col(ans)
  ans <- matrix(beta(k-1-(ii-1)+(jj-1), ii)/beta(k-1-(ii-1), ii), nrow = nrow(ans))
  ans
}

qda_params1 <- list(
  k = 100,
  r_train = 20,
  r_test = 1,
  p = 2,
  wishart_par = 3,
  snr = 20
)

qda_params2 <- modifyList(qda_params1, list(k = 1000))
t1 <- proc.time()
res <- do.call(simulate_qda, qda_params2)
mc_rates <- sapply(1:1000,
                   function(k) {
                     resample_misclassification(res$plikes, res$i_chosen, k)
                   })
proc.time() - t1
rm(res)

