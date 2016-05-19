#### Do simulations k = 2 to 50
source("extrapolation/simulation_source.R")
K <- 50
synth_data <- gen.data(p = 10, sigma = 1, K = K, r1 = 20, r2 = 30)
tabs <- list()
ks <- 3:50
for (k in ks) {
  tab <- build_extrapolation_table(synth_data, k, K)
  tabs <- c(tabs, list(tab))
}

