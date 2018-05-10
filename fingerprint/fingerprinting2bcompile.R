remove.sub332 <- logical()
remove.sub332["All_Sub_REST1_errts.mat"] <- TRUE
remove.sub332["All_Sub_REST2_errts.mat"] <- TRUE
remove.sub332["All_Sub_EMOTION_errts.mat"] <- TRUE
remove.sub332["All_Sub_GAMBLING_errts.mat"] <- FALSE
remove.sub332["All_Sub_LANGUAGE_errts.mat"] <- TRUE
remove.sub332["All_Sub_MOTOR_errts.mat"] <- FALSE
remove.sub332["All_Sub_RELATIONAL_errts.mat"] <- TRUE
remove.sub332["All_Sub_SOCIAL_errts.mat"] <- TRUE
remove.sub332["All_Sub_WM_errts.mat"] <- TRUE

tab <- matrix(1, 9, 9)
rownames(tab) <- sapply(sapply(sapply(names(remove.sub332), substr, start = 9, stop = 40), strsplit, split = "_"), `[[`, 1)
colnames(tab) <- rownames(tab)

tab2 <- tab

for (i in 1:8) {
  for (j in (i+1):9) {
    fname <- paste0("fingerprint/results/results", i, "_", j, ".txt")
    temp <- read.table(fname, header = FALSE)
    tab[i, j] <- temp[1,2]
    tab[j, i] <- temp[1,2]
    tab2[i, j] <- temp[2,2]
    tab2[j, i] <- temp[2,2]
    
  }
}


write.table(tab, file = "fingerprint/cor_ident.csv", sep = ",")
write.table(tab2, file = "fingerprint/KL_ident.csv", sep = ",")

pdf("fingerprint/cor_vs_KL.pdf", width = 15, height = 15)
plot(tab[row(tab) < col(tab)], tab2[row(tab) < col(tab)], xlab = "cor_ident", ylab = "KL_ident")
abline(0, 1)
for (i in 1:8) {
  for (j in (i+1):9) {
    text(tab[i,j], tab2[i, j] - 0.005, paste(rownames(tab)[i], "->", rownames(tab)[j]), cex = 0.6)
  }
}
dev.off()
