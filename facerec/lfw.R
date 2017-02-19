library(class)
library(pracma)

x <- read.csv("facerec/reps.csv", header = FALSE)
x <- as.matrix(x)
dim(x)
labs <- read.csv("facerec/labels.csv", header = FALSE, stringsAsFactors = FALSE)
dim(labs)
labs[, 2] <- sapply(labs[, 2], substring, 19)
y <- labs[, 1]
rownames(x) <- paste0(y, "_", (1:length(y)) %% 11)
save(x, y, file = "facerec/lfw.RData")
