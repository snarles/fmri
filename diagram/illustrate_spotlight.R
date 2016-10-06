library(png)
library(pracma)

## loading pics

icons <- list()
icons[[1]] <- readPNG("diagram/k_pic1.png")
icons[[2]] <- readPNG("diagram/k_pic2.png")
br <- readPNG("diagram/brain100.png")

## figure out that ellipse

plot(NA, NA, xlim = c(0, 1), ylim = c(0, 1))
rasterImage(br, 0, 0, 1, 1)
locator(1)

## ((col - 47)/28)^2 + ((row - 50)/38)^2

## generation of patterns

mtrx <- zeros(100, 100)
rc <- cbind(as.numeric(row(mtrx)), as.numeric(col(mtrx)))
dm <- distmat(rc, rc)

## drawing pics

cols <- c("red", "blue")
(alter <- rep(1:2, 5))

plot(NA, NA, xlim = c(-0.5, 10.5), ylim = c(0, 10), axes = FALSE, ann = FALSE)

for (i in 1:10) {
  j <- alter[i]
  polygon(c(0.25+i-1, 0.75+i-1, 0.75+i-1, 0.25+i-1), c(9.25, 9.25, 9.75, 9.75), 
          col = cols[j], border = cols[j])
  rasterImage(icons[[j]], 0.3+i-1, 9.3, 0.7+i-1, 9.7)
}

for (i in 1:10) {
  rasterImage(br, i-1+.1, 8, i-.1, 9)
}



