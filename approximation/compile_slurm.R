seq1 <- nits*10-9
seq2 <- nits * 10
str2 <- "/data/MLcore/temp/sim7"
lf <- paste0(str2, "_", seq1, "_", seq2, ".RData")
head(lf)

resZ <- list()
for (f in lf) {
  load(f)
  resZ <- c(resZ, res)
}
length(resZ)

res <- resZ

saveRDS(res, file = "approximation/raw_7_results.rds")
