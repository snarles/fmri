library(class)
library(pracma)
library(lineId)
source("extrapolation/ku_source.R")

# x <- read.csv("facerec/reps.csv", header = FALSE)
# x <- as.matrix(x)
# dim(x)
# labs <- read.csv("facerec/labels.csv", header = FALSE, stringsAsFactors = FALSE)
# dim(labs)
# labs[, 2] <- sapply(labs[, 2], substring, 19)
# y <- labs[, 1]
# rownames(x) <- labs[, 2]
# save(x, y, file = "facerec/lfw.RData")
load("facerec/lfw.RData")

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

set.seed(0)
nsub <- 400
clsub <- sample(ncol(split_ind_m), nsub)
tr_inds <- split_ind_m[1, clsub]
te_inds <- split_ind_m[2, clsub]
xtr <- xf[tr_inds, ]
xte <- xf[te_inds, ]
ytr <- yf[tr_inds]
yte <- yf[te_inds]

liks <- -pdist2(xte, xtr)
accs <- 1 - resample_misclassification(liks, 1:nsub, 1:nsub)

plot(accs, type = "l", ylim = c(0, 1), xlab = "faces", ylab = "accuracy",
     main = "Subset of 400", cex.lab = 1.5)

## full acc

tr_inds <- split_ind_m[1, ]
te_inds <- split_ind_m[2, ]
xtr <- xf[tr_inds, ]
xte <- xf[te_inds, ]
ytr <- yf[tr_inds]
yte <- yf[te_inds]

liks <- -pdist2(xte, xtr)
accs_full <- 1 - resample_misclassification(liks, 1:ncol(split_ind_m), 
                                       1:ncol(split_ind_m))


## fit spline model

nsplines <- 10000
acs <- 1 - accs
knts <- seq(0, 1, length.out = nsplines + 2)
knts <- knts[-c(1, nsplines + 2)]
MM <- spline1_moments(knts, 1:ncol(split_ind_m))
xmat <- spline1_moments(knts, 1:nsub)
bt <- nnls::nnls(xmat, acs)$x
xs <- seq(0, 1, 0.01)
plot(xs, spline1_dm(knts, xs) %*% bt, type = "l", main = "estimated D(u)", ylab = "D(u)", xlab = "u",
     cex.lab = 1.2)

acs_hat <- 1 - as.numeric(MM %*% bt)

pdf("facerec/acc_plot2.pdf")
plot(accs_full, type = "l", ylim = c(0, 1), xlab = "faces", ylab = "accuracy",
     main = "Full set (1672)", cex.lab = 1.5, col = "blue")
lines(1:nsub, accs, col = "black", lwd = 3)
lines(acs_hat, col = "red")
abline(v = nsub, lty = 2, col = "red")
legend(500, 0.3, c("sub400", "estimate", "full1672"), col = c("black", "red", "blue"), 
       lwd = c(2, 1, 1), cex = 1.5)
dev.off()

## multiple repeats

draw_sub_accs <- function(nsub) {
  clsub <- sample(ncol(split_ind_m), nsub)
  tr_inds <- split_ind_m[1, clsub]
  te_inds <- split_ind_m[2, clsub]
  xtr <- xf[tr_inds, ]
  xte <- xf[te_inds, ]
  ytr <- yf[tr_inds]
  yte <- yf[te_inds]
  liks <- -pdist2(xte, xtr)
  accs <- 1 - resample_misclassification(liks, 1:nsub, 1:nsub)
  nsplines <- 10000
  acs <- 1 - accs
  knts <- seq(0, 1, length.out = nsplines + 2)
  knts <- knts[-c(1, nsplines + 2)]
  MM <- spline1_moments(knts, 1:ncol(split_ind_m))
  xmat <- spline1_moments(knts, 1:nsub)
  bt <- nnls::nnls(xmat, acs)$x
  xs <- seq(0, 1, 0.01)
  acs_hat <- pmax(1 - as.numeric(MM %*% bt), 0)
  acs_hat
}

nreps <- 20
nsubs <- c(50, 100, 200, 400, 800)
#nsub <- 400

for (nsub in nsubs) {
  pdf(paste0("facerec/sub_", nsub, ".pdf"))
  plot(accs_full, type = "l", ylim = c(0, 1), xlab = "faces", ylab = "accuracy",
       main = paste(nsub), cex.lab = 1.5, col = "blue")
  abline(v = nsub, lty = 2)
  for (i in 1:nreps) lines(draw_sub_accs(nsub), col = "red")
  dev.off()
}

nreps <- 50
accs_final <- list()
for (nsub in nsubs) {
  accs <- sapply(1:nreps, function(i) min(draw_sub_accs(nsub)))
  accs_final[[paste(nsub)]] <- accs
}

boxplot(accs_final, ylim = c(0, 1), main = "Predicted accuracy (1672)")
abline(h = min(accs_full), col = "blue")

saveRDS(accs_final, file = "lfw_sub_preds.rds")
