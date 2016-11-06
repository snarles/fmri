source("info_theory_sims/sim4source.R")


Bmat <- 4 * eye(4)
compute_mi(Bmat)
run_simulation(Bmat, m.folds = 200, k.each = 10, r.each = 100, r.train = 50)





