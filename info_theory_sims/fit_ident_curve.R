## find the I_implied for classification curve

library(lineId)

fit_I_to_curve <- function(acs, ks = 1:length(acs), 
                           kmax = max(ks), nits = 10, i_max = 10,
                           wt_exp = 0) {
  wts <- ks^wt_exp
  i_lb <- 0
  i_ub <- i_max
  i_cur <- (i_lb + i_ub)/2
  for (i in 1:nits) {
    yh <- 1 - piK(sqrt(2 * i_cur), ks)
    dyh <- -d_piK(sqrt(2 * i_cur), ks)
    (of <- sum(wts * (acs - yh)^2))
    (deriv <- sum(wts * (yh - acs)*dyh))
    #print(c(i_cur = i_cur, of = of, deriv = deriv))
    if (deriv > 0) {
      ## decrease i_cur
      i_ub <- i_cur
      i_cur <- (i_lb + i_ub)/2
    } else {
      ## increase i_cur
      i_lb <- i_cur
      i_cur <- (i_lb + i_ub)/2
    }
  }
  i_cur
}

acs_curve <- function(i_implied, ks) {
  1 - piK(sqrt(2 * i_implied), ks)
}