subAcc <- function(pmat, cl) {
  cl_assigned <- apply(pmat, 2, function(v) which(v==max(v))[1])
  sum(cl == cl_assigned)/length(cl)
}

getUs <- function(pmat, ncl, ny) {
  rankconv <- (apply(pmat, 2, rank) - 0.5)/nrow(pmat)
  Us <- list()
  for (i in 1:ncl) {
    Us[[i]] <- rankconv[i, 1:ny + (ny) * (i-1)]  
  }
  Us <- do.call(c, Us)
  Us
}

getYs <- function(pmat, ncl, ny) {
  rankconv <- apply(pmat, 2, rank)
  Us <- list()
  for (i in 1:ncl) {
    Us[[i]] <- rankconv[i, 1:ny + (ny) * (i-1)]  
  }
  Us <- do.call(c, Us)
  Us
}


pmat <- read.table('/home/snarles/github/predict_test_error/lala.txt', header = FALSE)
saveRDS(pmat, "extrapolation/cifar100.rds")
cl <- rep(1:100, each = 100)

cl_assigned <- apply(pmat, 2, function(v) which(v==max(v))[1])
sum(cl == cl_assigned)/length(cl)

sa <- subAcc(pmat[1:30, 1:3000], cl[1:3000])
fa <- subAcc(pmat, cl)

Us <- getUs(pmat, 100, 100)

hist(getUs(pmat, 100, 100))

# subxs <- sort(sample(100, 30))
# subp <- pmat[subxs, cl %in% subxs]
# subus <- getUs(subp, 30, 100)

library(lineId)
(ihat <- Ihat_LI(1 - mean(sa), 30))
(fa_i <- 1 - lineId::piK(sqrt(2 * ihat), 100))

### Try MLE/constrained MLE
sub_us <- getYs(pmat[1:30, 1:3000], ncl = 30, ny = 100)
mle_est <-  res_mixtools(sub_us, 30)
pseq <- seq(0, 1, 1/300)
k <- 30
cm <- cons_mle_est(Ys, k, pseq, 0.01)
(fa_m <- est_moment(mle_est, 100))
(fa_c <- cm_est_moment(cm, 100))

c(fa, fa_i, fa_m, fa_c)
