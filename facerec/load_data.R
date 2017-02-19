library(class)
library(pracma)

x <- read.csv("att_faces_aligned/reps.csv", header = FALSE)
x <- as.matrix(x)
View(x)
dim(x)
labs <- read.csv("att_faces_aligned/labels.csv", header = FALSE, stringsAsFactors = FALSE)
dim(labs)
help(substr)
labs[, 2] <- sapply(labs[, 2], substring, 19)
y <- labs[, 1]
rownames(x) <- paste0(y, "_", (1:length(y)) %% 11)

floor(pdist(x[y <= 2, ]) * 10)

## distance matrix
image(pdist(x[y <= 10, ]))


## built test set
test_inds <- sapply(1:40, function(x) sample(which(labs[, 1]==x), 1))
x.te <- x[test_inds, ]
x.tr <- x[-test_inds, ]
y.te <- y[test_inds]
y.tr <- y[-test_inds]

## alternatively, build training set
tr_inds <- sapply(1:40, function(x) sample(which(labs[, 1]==x), 1, replace = FALSE))
x.te <- x[-tr_inds, ]
x.tr <- x[tr_inds, ]
y.te <- y[-tr_inds]
y.tr <- y[tr_inds]


## success of knn???
knn_errs <- numeric()
k <- 1
#for (k in 1:5) {
  res <- knn(train = x.tr, test = x.te, cl = y.tr, k = k)
  knn_errs[paste0("k_", k)] <- mean(res!=y.te)
#}
knn_errs

