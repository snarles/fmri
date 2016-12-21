## find the I_implied for classification curve

library(lineId)

fit_I_to_curve <- function(acs, ks = 1:length(acs), 
                           kmax = max(ks), nits = 3, i_max = 10,
                           wts = sqrt(kmax/ks)) {
  i_lb <- 0
  i_ub <- i_max
  i_cur <- (i_lb + i_ub)/2
  yh <- 1 - sapply(ks, function(k) piK(sqrt(2 * i_cur), k))
  
}