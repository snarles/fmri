## find the I_implied for classification curve

library(lineId)

fit_I_to_curve <- function(acs, ks = 1:length(acs), 
                           kmax = max(ks), nits = 3, i_max = 10,
                           wts = sqrt(kmax/ks)) {
  i_lb <- 0
  i_ub <- i_max
  i_cur <- (i_lb + i_ub)/2
  yh <- 1 - piK(sqrt(2 * i_cur), ks)
  dyh <- -d_piK(sqrt(2 * i_cur), ks)
}


numDeriv::grad(function(v) piK(v, 3), 0.5)
d_piK(0.5, 3)

numDeriv::grad(function(v) piK(v, 4), 0.5)
d_piK(0.5, 3:5)
