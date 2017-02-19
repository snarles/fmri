library(class)
library(pracma)
library(lineId)
source("extrapolation/ku_source.R")

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

nsub <- 1600
clsub <- sample(ncol(split_ind_m), nsub)
tr_inds <- split_ind_m[1, clsub]
te_inds <- split_ind_m[2, clsub]
xtr <- xf[tr_inds, ]
xte <- xf[te_inds, ]
ytr <- yf[tr_inds]
yte <- yf[te_inds]

liks <- -pdist2(xte, xtr)
accs <- 1 - resample_misclassification(liks, 1:nsub, 1:nsub)

plot(accs, type = "l", ylim = c(0, 1), xlab = "faces", ylab = "accuracy")
acs_hat <- 1 - spline_extrap(1 - accs, 200, nsplines = 10000)
lines(acs_hat, col = "red")

res <- knn(train = xtr, test = xte, cl = ytr, k = 1)
(knn_err <- mean(res!=yte))




