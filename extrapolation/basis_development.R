## develop a good basis to use

source("extrapolation/ku_source.R")



load("approximation/sim_large5_k5_raw.RData", verbose = TRUE)

plot(kref, accs_subs[5, ], type = "l")
plot(Ktarg, accsZ[5, ], type = "l")
