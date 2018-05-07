str1 <- "Rscript /data/MLcore/fmri/approximation/mc_gaussian_sim7_batch.R" 
str2 <- "/data/MLcore/temp/sim7"
nits <- 1:1000
sink("approximation/swarm_sim7.swarm")
cat(paste(paste(str1, nits*10-9, nits * 10, str2), collapse = "\n"))
sink()