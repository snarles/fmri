## linear simulation


library(pracma)
source("idloss/mi_source.R")
source("idloss/idLoss.R")
n <- 1000
p <- 2
sgma <- 0.1

RES <- list()
#sink("idloss/temp.txt")
#d_ex <- 0
for (d_ex in 100 * 1:6) {

  set.seed(0)
  
  X0 <- randn(n, p)
  Y0 <- X0 + sgma * randn(n, p)
  ## true MI
  (mi_true <- p/2 * log(1 + sgma^(-2)))
  
  bmat1 <- 0.1 * randn(p) + eye(p)
  bmat2 <- 0.1 * randn(p) + eye(p)
  X <- X0 %*% bmat1
  Y <- Y0 %*% bmat2
  colnames(X) <- paste0("X", 1:p)
  colnames(Y) <- paste0("Y", 1:p)
  pairs(cbind(X, Y))
  
  if (d_ex > 0) {
    X <- cbind(X, randn(n, d_ex))
    #Y <- cbind(Y, randn(n, d_ex))
  }
  colnames(X) <- paste0("X", 1:ncol(X))
  
  
  
  ## nn estimate
  (mi_nn <- nn_mi(X, Y)[1])
  
  ## id loss using OLS
  k <- 10
  (idl <- id_cv_loss(X, Y, k, mc.reps = 1000))
  (mi_li <- lineId::Ihat_LI(idl, k))
  (mi_np_ols_10 <- lineId::aba_to_mi_lower(k, 1 - idl))
  
  k <- 20
  (idl <- id_cv_loss(X, Y, k, mc.reps = 1000))
  (mi_li <- lineId::Ihat_LI(idl, k))
  (mi_np_ols_20 <- lineId::aba_to_mi_lower(k, 1 - idl))
  
  ## id loss using enet
  k <- 10
  (idl <- id_cv_loss(X, Y, k, mc.reps = 1000, fitter = fitter_enet, alpha = 0.5))
  (mi_li <- lineId::Ihat_LI(idl, k))
  (mi_np_enet_10 <- lineId::aba_to_mi_lower(k, 1 - idl))
  
  k <- 20
  (idl <- id_cv_loss(X, Y, k, mc.reps = 1000, fitter = fitter_enet, alpha = 0.5))
  (mi_li <- lineId::Ihat_LI(idl, k))
  (mi_np_enet_20 <- lineId::aba_to_mi_lower(k, 1 - idl))
  
  
  
  ## id loss using RF
  # k <- 10
  # (idl <- id_cv_loss(X, Y, k, mc.reps = 50, fitter = fitter_rf))
  # (mi_li <- lineId::Ihat_LI(idl, k))
  # (mi_np_rf_10 <- lineId::aba_to_mi_lower(k, 1 - idl))
  # 
  # k <- 20
  # (idl <- id_cv_loss(X, Y, k, mc.reps = 50, fitter = fitter_rf))
  # (mi_np_rf_20 <- lineId::aba_to_mi_lower(k, 1 - idl))
  # 
  reso <- c(mi_true = mi_true, mi_nn = mi_nn, 
       mi_np_ols_10 = mi_np_ols_10, mi_np_ols_20 = mi_np_ols_20,
       mi_np_enet_10 = mi_np_enet_10, mi_np_enet_20 = mi_np_enet_20)
  #,mi_np_rf_10 = mi_np_rf_10, mi_np_rf_20 = mi_np_rf_20)
  RES[[paste0("d", d_ex)]] <- reso
}
res2 <- data.frame(d_ex = 1:6 * 100, do.call(rbind, RES))
res2
saveRDS(res2, "idloss/sim2a_results.rds")


## MAKE PLOTS

res1 <- readRDS("idloss/sim2_res.rds")
res2 <- readRDS("idloss/sim2_res2.rds")
res3 <- readRDS("idloss/sim2a_results.rds")

View(res1)

ress <- res1
plot(NA, NA, ylim = c(0, 5), xlim = c(0, 22), xlab = "dimension of X",
     ylab = "estimated MI")
abline(h = res1$mi_true[1], lwd = 2, lty = 2)
cols <- c(nn = rgb(0.8,0,0), ols10 = rgb(0,0,0.4),
          ols20 = rgb(0,0,0.6), enet10 = rgb(0,0.4,0),
          enet20 = rgb(0,0.6,0))
lines(ress$d_ex + 2, ress$mi_nn, col = cols["nn"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_ols_10, col = cols["ols10"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_ols_20, col = cols["ols20"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_enet_10, col = cols["enet10"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_enet_20, col = cols["enet20"], lwd = 2, type = "o")
legend(13, 4, lwd = 2, 
       col = c("black", cols), legend = c("true", "nonpar", "OLS (k=10)", "OLS (k=20)",
                              "enet (k=10)", "enet (k=20)"),
       lty = c(2, 1,1,1,1,1))


plot(NA, NA, ylim = c(4.1, 4.7), xlim = c(0, 22), xlab = "dimension",
     ylab = "estimated MI")
abline(h = res1$mi_true[1], lwd = 2)
cols <- c(nn = grey(0.2), ols10 = rgb(0,0,0.4),
          ols20 = rgb(0,0,0.6), enet10 = rgb(0,0.4,0),
          enet20 = rgb(0,0.6,0))
lines(res1$d_ex + 2, res1$mi_nn, col = cols["nn"], lwd = 2, type = "o")
lines(res1$d_ex + 2, res1$mi_np_ols_10, col = cols["ols10"], lwd = 2, type = "o")
lines(res1$d_ex + 2, res1$mi_np_ols_20, col = cols["ols20"], lwd = 2, type = "o")
lines(res1$d_ex + 2, res1$mi_np_enet_10, col = cols["enet10"], lwd = 2, type = "o")
lines(res1$d_ex + 2, res1$mi_np_enet_20, col = cols["enet20"], lwd = 2, type = "o")
legend(5, 4.3, lwd = 2, 
       col = cols[1:3], legend = c("nn", "ols10", "ols20"))
legend(15, 4.25, lwd = 2, 
       col = cols[4:5], legend = c("enet10", "enet20"))



res2$d_ex
ress <- res2
plot(NA, NA, ylim = c(4, 4.8), xlim = c(0, 202), xlab = "dimension of X",
     ylab = "estimated MI")
abline(h = res1$mi_true[1], lwd = 2)
cols <- c(nn = grey(0.2), ols10 = rgb(0,0,0.4),
          ols20 = rgb(0,0,0.6), enet10 = rgb(0,0.4,0),
          enet20 = rgb(0,0.6,0))
#lines(ress$d_ex + 2, ress$mi_nn, col = cols["nn"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_ols_10, col = cols["ols10"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_ols_20, col = cols["ols20"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_enet_10, col = cols["enet10"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_enet_20, col = cols["enet20"], lwd = 2, type = "o")
legend(50, 4.22, lwd = 2, 
       col = cols[2:3], legend = c("ols10", "ols20"))
legend(125, 4.22, lwd = 2, 
       col = cols[4:5], legend = c("enet10", "enet20"))


res3$d_ex
ress <- res3
plot(NA, NA, ylim = c(3, 4.8), xlim = c(0, 602), xlab = "dimension of X",
     ylab = "estimated MI")
abline(h = res1$mi_true[1], lwd = 2, lty = 2)
cols <- c(nn = grey(0.2), ols10 = rgb(0,0,0.4),
          ols20 = rgb(0,0,0.6), enet10 = rgb(0,0.4,0),
          enet20 = rgb(0,0.6,0))
#lines(ress$d_ex + 2, ress$mi_nn, col = cols["nn"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_ols_10, col = cols["ols10"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_ols_20, col = cols["ols20"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_enet_10, col = cols["enet10"], lwd = 2, type = "o")
lines(ress$d_ex + 2, ress$mi_np_enet_20, col = cols["enet20"], lwd = 2, type = "o")
legend(50, 3.75, lwd = 2, 
       col = c("black", cols[-1]), legend = c("true", "OLS (k=10)", "OLS (k=20)",
                                          "enet (k=10)", "enet (k=20)"),
       lty = c(2, 1,1,1,1))

