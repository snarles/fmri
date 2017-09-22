str1 <- "Rscript /data/MLcore/fmri/approximation/mc_gaussian_sim6a_batch.R" 
str2 <- "/data/MLcore/temp/sim6a"
nits <- 1:400
sink("approximation/swarm_sim6a.swarm")
cat(paste(paste(str1, nits*10-9, nits * 10, str2), collapse = "\n"))
sink()