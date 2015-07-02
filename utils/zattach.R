zattach <- function(ll) {
  for (i in 1:length(ll)) {
    assign(names(ll)[i], ll[[i]], envir=globalenv())
  }
}

#dots <- function(...) {
#  eval(substitute(alist(...)))
#}

do.call2 <- function(f, ll, ...) {
  dots <- eval(substitute(alist(...)))
  ll <- modifyList(ll, dots)
  do.call(f, ll)
}