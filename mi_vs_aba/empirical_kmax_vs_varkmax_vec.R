library(AlgDesign)


compute_moments <- function(d, alpha, k, mc.reps = 1e4) {
  maxstat <- sapply(1:mc.reps, function(i) {
    aa <- matrix(rgamma(d * k, alpha), k, d)
    rs <- rowSums(aa)
    aa[rs == 0, ] <- 1
    rs <- rowSums(aa)
    aa <- aa/rowSums(aa)
    lala <- apply(aa, 2, max)/k
    sum(lala)
  })
  c(mean(maxstat), sd(maxstat))
}

####
## k = 3
####

## d = 2
compute_moments(2, 2, 3)    # 0.4617648 0.0609609
compute_moments(2, 8e-3, 3) # 0.5835685 0.1420055
compute_moments(2, 4e-3, 3) # 0.5837331 0.1428442
compute_moments(2, 2e-3, 3) # 0.5858483 0.1344736


## d = 3
compute_moments(3, 2, 3)    # 0.4894503 0.0593871
compute_moments(3, 0.1, 3)  # 0.6865047 0.1482515
compute_moments(3, 8e-3, 3) # 0.7024501 0.1858227
compute_moments(3, 4e-3, 3) # 0.7025812 0.1868771
compute_moments(3, 2e-3, 3) # 0.7054425 0.1810342

## d = 4
compute_moments(4, 5e-3, 3)

