library(class)
library(pracma)

x <- read.csv("facerec/reps.csv", header = FALSE)
x <- as.matrix(x)
dim(x)
labs <- read.csv("facerec/labels.csv", header = FALSE, stringsAsFactors = FALSE)
dim(labs)
labs[, 2] <- sapply(labs[, 2], substring, 19)
y <- labs[, 1]
rownames(x) <- labs[, 2]
save(x, y, file = "facerec/lfw.RData")

## get 2-repeat data
yc <- table(y)
y2s <- as.numeric(names(yc)[yc > 1])
head(y2s)
filt <- y %in% y2s
xf <- x[filt, ]; yf <- y[filt]

split_inds <- sapply(unique(yf), function(x) sample(which(yf==x), 2))
length(split_inds)
split_ind_m <- matrix(split_inds, nrow = 2)

dim(split_ind_m)

nsub <- 100
clsub <- sample(ncol(split_ind_m), nsub)
tr_inds <- split_ind_m[1, clsub]
te_inds <- split_ind_m[2, clsub]
xtr <- xf[tr_inds, ]
xte <- xf[te_inds, ]
ytr <- yf[tr_inds]
yte <- yf[te_inds]

res <- knn(train = xtr, test = xte, cl = ytr, k = 1)
(knn_err <- mean(res!=yte))



