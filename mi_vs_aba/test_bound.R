source("mi_vs_aba/alpha_vs_beta.R")
source("mi_vs_aba/code1.R")

library(pracma)
library(lineId)
p <- 200
mat <- exp(randn(p)) + 1000 * eye(p)
mat <- colnorm(mat, nits = 20)
#rowSums(mat)

#mi_mat(mat)
(i_true <- mi_mat2(mat))
k <- 6
(aba <- aba_mat(mat, k))
(abas <- aba_mat_trials(mat, k, ntr = 20))
#aba_mat_naive(mat, k = 2, nits = 1e4)
(i_np <- find_par_aba(aba, k)$res[1])
find_par_aba(min(abas), k)
find_par_aba(max(abas), k)
(i_li <- lineId::Ihat_LI(1-aba, k))
cbind(i_true, i_np, i_li)
