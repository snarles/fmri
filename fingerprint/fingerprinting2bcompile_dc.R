fprefix <- "/data/HCP_preproc/fingerprinting/Results/DC_Matrices/"
nms <- list.files(fprefix)


task_names <- sapply(sapply(sapply(nms, substr, start = 9, stop = 40), strsplit, split = "_"), `[[`, 1)
gsrs <- grep(pattern = 'GSR', x = nms)
not_gsrs <- grep(pattern = 'GSR', x = nms, invert = TRUE)

tab <- matrix(1, 18, 18)
#rownames(tab) <- sapply(sapply(sapply(nms, substr, start = 9, stop = 40), strsplit, split = "_"), `[[`, 1)
rownames(tab) <- nms
colnames(tab) <- rownames(tab)

tab2 <- tab

for (i in 1:17) {
  for (j in (i+1):18) {
    fname <- paste0("fingerprint/results/results_dc", i, "_", j, ".txt")
    temp <- read.table(fname, header = FALSE)
    tab[i, j] <- temp[1,2]
    tab[j, i] <- temp[3,2]
    tab2[i, j] <- temp[2,2]
    tab2[j, i] <- temp[4,2]
    
  }
}


write.table(tab, file = "fingerprint/cor_ident_dc.csv", sep = ",")
write.table(tab2, file = "fingerprint/KL_ident_dc.csv", sep = ",")

tab_gsr <- tab[gsrs, gsrs]
tab2_gsr <- tab2[gsrs, gsrs]
rownames(tab_gsr) <- task_names[gsrs]
colnames(tab_gsr) <- task_names[gsrs]
rownames(tab2_gsr) <- task_names[gsrs]
colnames(tab2_gsr) <- task_names[gsrs]

tab_not_gsr <- tab[not_gsrs, not_gsrs]
tab2_not_gsr <- tab2[not_gsrs, not_gsrs]
rownames(tab_not_gsr) <- task_names[not_gsrs]
colnames(tab_not_gsr) <- task_names[not_gsrs]
rownames(tab2_not_gsr) <- task_names[not_gsrs]
colnames(tab2_not_gsr) <- task_names[not_gsrs]

tab <- tab_gsr; tab2 <- tab2_gsr
pdf("fingerprint/cor_vs_KL_dc_gsr.pdf", width = 15, height = 15)
plot(tab[row(tab) != col(tab)], tab2[row(tab) != col(tab)], xlab = "cor_ident", ylab = "KL_ident")
abline(0, 1)
for (i in 1:9) {
  for (j in 1:9) {
    if (i != j)
      text(tab[i,j], tab2[i, j] - 0.005, paste(rownames(tab)[i], "->", rownames(tab)[j]), cex = 0.6)
  }
}
dev.off()

tab <- tab_not_gsr; tab2 <- tab2_not_gsr
pdf("fingerprint/cor_vs_KL_dc_not_gsr.pdf", width = 15, height = 15)
plot(tab[row(tab) != col(tab)], tab2[row(tab) != col(tab)], xlab = "cor_ident", ylab = "KL_ident")
abline(0, 1)
for (i in 1:9) {
  for (j in 1:9) {
    if (i != j)
      text(tab[i,j], tab2[i, j] - 0.005, paste(rownames(tab)[i], "->", rownames(tab)[j]), cex = 0.6)
  }
}
dev.off()

tab <- tab_not_gsr; tab2 <- tab_gsr
pdf("fingerprint/cor_not_gsr_vs_gsr.pdf", width = 15, height = 15)
plot(tab[row(tab) != col(tab)], tab2[row(tab) != col(tab)], xlab = "notGSR", ylab = "GSR")
abline(0, 1)
for (i in 1:9) {
  for (j in 1:9) {
    if (i != j)
      text(tab[i,j], tab2[i, j] - 0.005, paste(rownames(tab)[i], "->", rownames(tab)[j]), cex = 0.6)
  }
}
dev.off()

tab <- tab2_not_gsr; tab2 <- tab2_gsr
pdf("fingerprint/KL_not_gsr_vs_gsr.pdf", width = 15, height = 15)
plot(tab[row(tab) != col(tab)], tab2[row(tab) != col(tab)], xlab = "notGSR", ylab = "GSR")
abline(0, 1)
for (i in 1:9) {
  for (j in 1:9) {
    if (i != j)
      text(tab[i,j], tab2[i, j] - 0.005, paste(rownames(tab)[i], "->", rownames(tab)[j]), cex = 0.6)
  }
}
dev.off()

