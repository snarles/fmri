## birthday collision probs

library(AlgDesign)

## naive count
naive_count <- function(k, d) {
  mat <- gen.factorial(rep(d, k), center = FALSE)
  ct <- apply(mat, 1, function(v) length(unique(v)))
  table(ct)
}

## theoretical table
build_table <- function(k, d) {
  maxl <- pmin(k, d)
  ctab <- matrix(0, d, maxl)
  ftab <- ctab
  
  for (dd in 1:d) {
    for (ll in 1:pmin(k, dd)) {
      if (ll == 1) {
        ctab[dd, 1] <- dd
        ftab[dd, 1] <- ctab[dd, 1]
      } else {
        ctab[dd, ll] <- choose(dd, ll) * (ll^k - ftab[ll, ll-1])
        ftab[dd, ll] <- ftab[dd, ll-1] + ctab[dd, ll]
      }
    }
  }

  list(ans = ctab[nrow(ctab), ], ctab = ctab, ftab = ftab)
}


moments_k <- function(k, dmax, naive = FALSE, mc.reps = 1e4) {
  if (naive) {
    ans <- matrix(0, dmax, 2)
    for (d in 1:dmax) {
      samps <- sapply(1:mc.reps, function(i) {
        v <- sample(d, k, TRUE)
        length(unique(v))
      })
      ans[d, ] <- c(mean(samps/k), sd(samps/k))
    }
    return(ans)
  }
  tab <- build_table(k, dmax)$ctab
  ans <- matrix(0, dmax, 2)
  for (i in 1:dmax) {
    v <- tab[i, ]
    v <- v/sum(v)
    vals <- (1:length(v))/k
    ev <- sum(vals * v)
    ev2 <- sum(vals^2 * v)
    va <- ev2 - (ev)^2
    ans[i, ] <- c(ev, sqrt(va))
  }
  ans
}

# naive_count(3, 4)
# build_table(3, 4)$ans
# 
# 
# naive_count(5, 4)
# build_table(5, 4)$ans
# 
# 
# naive_count(5, 5)
# build_table(5, 5)$ans
# 
# naive_count(5, 6)
# build_table(5, 6)$ans
# 
# naive_count(6, 6)
# build_table(6, 6)$ans
# 
# build_table(6, 8)$ans
# build_table(6, 10)$ctab[8, ]



# moments_k(2, 3)
# moments_k(2, 3, TRUE)
# 
# moments_k(3, 6)
# moments_k(3, 6, TRUE)
# 
# moments_k(4, 8)
# moments_k(4, 8, TRUE)
# 
# moments_k(5, 10)
