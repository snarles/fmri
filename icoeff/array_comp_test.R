## faster way to get colsums of array??

ad <- c(1000, 400, 300)
arr1 <- array(rnorm(prod(ad)), dim = ad)

t1 <- proc.time()
res1 <- apply(arr1, c(2, 3), sum)
proc.time() - t1
dim(res1)

t2 <- proc.time()
mat <- matrix(arr1, ad[1], ad[2] * ad[3])
res2 <- matrix(colSums(mat), ad[2], ad[3])
proc.time() - t2

sum(abs(res1 - res2))
