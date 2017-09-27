str1 <- "Rscript /data/MLcore/fmri/approximation/mc_gaussian_sim5_collect_data.R" 
str2 <- "/data/MLcore/temp/sim5cd"
nits <- 1:1000
sink("approximation/swarm_sim5cd.swarm")
cat(paste(paste(str1, nits*10-9, nits * 10, str2), collapse = "\n"))
sink()