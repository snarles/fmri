p2n <- function(x) as.numeric(sub("%","", x))/100

err1 <- read.table("rakesh/error_rates/err_rates.delim", sep = " ", header = TRUE)
View(err1)
err1$TrainErr <- err1$TrainErr/100
err1$TestErr <- err1$TestErr/100

err_knn400 <- read.table("rakesh/error_rates/knn400.txt", header = FALSE)
colnames(err_knn400) <- c("k", "tr", "te")
err_knn400[, 2] <- p2n(err_knn400[, 2])
err_knn400[, 3] <- p2n(err_knn400[, 3])

View(err_knn400)

lprobs <- list()
fl <- list.files("rakesh/error_rates/nnet_probs")
for (ff in fl) {
  tab <- read.table(paste0("rakesh/error_rates/nnet_probs/", ff), header = FALSE)
  lprobs[[ff]] <- tab
}

knnprobs <- list()
(fl <- list.files("rakesh/error_rates/knn_probs"))
for (ff in fl) {
  tab <- read.table(paste0("rakesh/error_rates/knn_probs/", ff), header = FALSE)
  knnprobs[[ff]] <- tab
}

save(err_knn400, err1, knnprobs, lprobs, file = "rakesh/converted1.rda")


####
##  New 100-class
####

lprobs <- list()
fl <- list.files("rakesh/sub_sub_runs/")
for (ff in fl) {
  tab <- read.table(paste0("rakesh/sub_sub_runs/", ff), header = FALSE)
  lprobs[[ff]] <- tab
}
saveRDS(lprobs, file = "rakesh/converted2.rds")

###
## Train_tel
###

errates <- read.table("rakesh/train_tel/error_rates", sep = " ", header = TRUE)
errates[, 3:4] <- 1 - errates[, 3:4]/100
colnames(errates)[3:4] <- c("tr_acc", "te_acc")
true_accs <- errates

lprobs <- list()
fl <- list.files("rakesh/train_tel/sub_sub_runs/")
for (ff in fl) {
  tab <- t(as.matrix(read.table(paste0("rakesh/train_tel/sub_sub_runs/", ff), header = FALSE)))
  lprobs[[ff]] <- tab
}

accs400 <- true_accs[true_accs$num_sub_classes==400, 4]
names(accs400) <- true_accs[true_accs$num_sub_classes==400, 1]
accs400

accs100 <- true_accs[true_accs$num_sub_classes==100, 4]
names(accs100) <- true_accs[true_accs$num_sub_classes==100, 1]
accs100


lp20 <- lprobs[c(1,3,6)]
lp100 <- lprobs[c(2,4,5)]
methods <- names(gt_accs)
names(lp100) <- methods
names(lp20) <- methods

save(true_accs, lprobs, lp20, lp100, methods, accs100, accs400, file = "rakesh/converted3.rda")
