ff <- scan("info_theory_paper/extrapolation_simple/temp.txt", character())
ff2 <- scan("info_theory_paper/extrapolation_simple/tempp.txt", character())

cmds <- sapply(1:length(ff), function(i)
  paste0("convert -units PixelsPerInch ",ff[i]," -density 72 ",ff2[i],"\n"))
sink("info_theory_paper/extrapolation_simple/temp2.txt")
for (c in cmds) cat(c)
sink()
