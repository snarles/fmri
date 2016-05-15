source("info_theory_sims/sim3source.R")

####
##  Identity case
####

allresults <- list()

## parallelization
mc.reps <- 1e4
mc.abe <- 1e2
abe.each <- 1e2
mcc <- 39
data.reps <- 750

## problem params
ss <- sqrt(seq(0, 200, by = 5))
ress <- array(0, dim = c(data.reps, 10, length(ss)))
p <- 10
mi_trues <- sapply(ss, function(s) mi_ident_case(p, s/sqrt(p), 1e5))
for (ii in 1:length(ss)) {
  mult <- ss[ii]/sqrt(p)
  Bmat <- mult * eye(p)
  (mi_true <- mi_ident_case(p, mult, 1e5))
  ## bayes LS
  k.each <- 20
  t1 <- proc.time()
  (est_ls <- get_abe(Bmat, k.each, abe.each, mc.abe, mcc))
  proc.time() - t1
  ## data params
  m.folds <- 1
  r.each <- 1000
  r.train <- floor(0.5 * r.each)
  (N = m.folds * k.each * r.each)
  # t1 <- proc.time()
  # run_simulation(Bmat, m.folds, k.each, r.each, r.train)
  # proc.time() - t1
  ## full-scale
  t1 <- proc.time()
  res <- run_simulations(Bmat, m.folds, k.each, r.each, r.train, mcc, data.reps)
  proc.time() - t1
  ress[, , ii] <- res
}
ress[, , 1]

li <- floor(0.1 * data.reps + 1); ui <- floor(0.9 * data.reps + 1)
m_i <- floor(0.5 * data.reps)
lowers <- t(apply(ress[,c(1:4, 8),], c(2, 3), function(v) sort(v)[li]))
uppers <- t(apply(ress[,c(1:4, 8),], c(2, 3), function(v) sort(v)[ui]))
meds <- t(apply(ress, c(2, 3), function(v) sort(v)[m_i]))
colnames(lowers) <- colnames(res)
colnames(uppers) <- colnames(res)

plot(NA, NA, xlab = "I", ylab = expression(hat(I)), xlim = c(0, max(mi_trues)),
     ylim = c(0, max(uppers)))
cols <- c(hsv(h = 0:3/4, s = 0.9, v = 0.5), hsv(0, 0, 0))
for (i in 1:5) {
  #lines(mi_trues, meds[, i], col = cols[i], lwd = 4)
  polygon(c(mi_trues, rev(mi_trues)), c(lowers[, i], rev(uppers[, i])), col = cols[i],
          border = cols[i])  
}
for (i in 1:5) {
  polygon(c(mi_trues, rev(mi_trues)), c(lowers[, i], rev(uppers[, i])), col = NA,
          border = cols[i], lwd = 3)  
}
abline(0, 1, lwd = 3, col = "white", lty = 2)

## save results
packet <- list(ss=ss, m.folds = m.folds,
               k.each = k.each, r.each = r.each, r.train = r.train,
               mi_true, est_ls = est_ls, ress = ress,
               mc.reps = mc.reps, mc.abe = mc.abe)
allresults <- c(allresults, list(packet))


save(allresults, file = 'info_theory_sims/fig4b.Rdata')