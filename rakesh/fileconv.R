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
