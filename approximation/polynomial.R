
check_divis <- function(x, n) {
  max(abs(x * n - round(x * n)))
}
display_div <- function(x, n) {
  nms <- names(x)
  x <- round(x * n)
  stx <- paste0(round(x * n), "/", n)
  names(stx) <- nms
  stx[x != 0]
}

mat_pows <- function(X, pows) {
  ans <- rep(1, dim(X)[1])
  for (i in 1:length(pows)) {
    ans <- ans * X[, i]^pows[i]
  }
  ans
}

## form matrix of all interactions to degree deg
form_poly_matrix <- function(X, deg) {
  p <- dim(X)[2]
  degs <- AlgDesign::gen.factorial(p + 1, deg)
  degs <- (degs - min(degs))/2
  tots <- t(apply(degs, 1, function(v) sapply(1:p, function(i) sum(v==i))))
  temp <- apply(tots, 1, function(v) paste(v, collapse = ","))
  filt <- match(sort(unique(temp)), temp)
  tots <- tots[filt, ]
  nn <- dim(tots)[1]
  res <- matrix(0, dim(X)[1], nn)
  nms <- character(nn)
  for (i in 1:nn) {
    tt <- tots[i, ]
    nms[i] <- paste(colnames(X)[tt > 0], tt[tt > 0], collapse = "", sep = "")
    res[, i] <- mat_pows(X, tt)
  }
  nms[1] <- "intercept"
  colnames(res) <- nms
  res
}

write_poly_formula <- function(X, coef) {
  p <- dim(X)[2]
  degs <- AlgDesign::gen.factorial(p + 1, deg)
  degs <- (degs - min(degs))/2
  tots <- t(apply(degs, 1, function(v) sapply(1:p, function(i) sum(v==i))))
  temp <- apply(tots, 1, function(v) paste(v, collapse = ","))
  filt <- match(sort(unique(temp)), temp)
  tots <- tots[filt, ]
  inds <- which(coef != 0)
  str <- character(length(inds))
  for (i in 1:length(inds)) {
    v <- coef[inds[i]]
    tt <- tots[inds[i], ]
    if (sum(tt) == 0) {
      str[i] <- paste(v)
    } else {
      polystrs <- colnames(X)[order(-tt)]
      pows <- tt[order(-tt)]
      polystrs <- polystrs[pows > 0]
      pows <- pows[pows > 0]
      polystrs[pows > 1] <- paste0(polystrs[pows > 1], "^", pows[pows > 1])
      str[i] <- paste(c(paste(v), polystrs), collapse = "*")
    }
  }
  paste0(str, collapse = " + ")  
}

