library(parallel)

avger <- function(f, mc.reps, mc.cores, avg, ...) {
  avgs <- mclapply(1:mc.reps, function(i) {
    set.seed(i)
    f(...)
  }, mc.cores =mc.cores)
  if (avg) return(mean(unlist(avgs)))
  avgs
}
