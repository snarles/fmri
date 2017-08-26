



generate_mat <- function(ny, nt, shape1 = 1) {
  mat <- matrix(rgamma(ny * nt, shape1), ny, nt)
  mat <- mat/rowSums(mat)
  mat
}

compute_moments <- function(mat, k, mc.reps = 100) {
  maxstat <- sapply(1:mc.reps, function(i) {
    aa <- mat[sample(nrow(mat), k, replace = TRUE), ]
    lala <- apply(aa, 2, max)/k
    sum(lala)
  })
  c(mean(maxstat), sd(maxstat))
}


ny <- 100
nt <- 100
shape1 <- 0.1
mat <- generate_mat(ny, nt, shape1)
rowSums(mat)
# image(mat)


mats <- lapply(1:1000, function(i) {
  generate_mat(ny, nt, exp(3 * rnorm(1)))
})

allmoms <- list()
for (k in 2:20) {
  moms <- do.call(rbind, lapply(mats, compute_moments, k = k))
  allmoms <- c(allmoms, list(moms))
  plot(moms, pch = ".")
  title(k)
}

plot(do.call(rbind, allmoms), pch = ".", col = "white")
for (i in 1:length(allmoms)) {
  points(allmoms[[i]], col = rainbow(length(allmoms))[i], pch = ".", cex = i/5)
}

maxes <- sapply(allmoms, function(j) max(j[, 2], na.rm = TRUE))
plot(maxes)
