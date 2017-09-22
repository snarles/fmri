seq1 <- nits*10-9
seq2 <- nits * 10
str2 <- "/data/MLcore/temp/sim6a"
lf <- paste0(str2, "_", seq1, "_", seq2, ".RData")

resZ <- list()
for (f in lf) {
  load(f)
  resZ <- c(resZ, res)
}
length(resZ)

res <- resZ
