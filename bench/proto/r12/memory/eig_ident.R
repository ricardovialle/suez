# Rscript eig_ident.R <lib> <tag>: kernel_contrast_eigen on test kernels -> out/eig_<tag>.rds
a = commandArgs(TRUE); suppressPackageStartupMessages(library(suez, lib.loc = a[1])); if (length(a) >= 3) options(suez.eigen_solver = a[3])
set.seed(3)
ks = list()
X = matrix(rnorm(300 * 40), 300); ks$rand = tcrossprod(X)                       # rank 40 -> many ties at 0
id = rep(1:150, each = 2)[sample(300)]; ks$zz = outer(id, id, "==") * 1          # default kernel
G = matrix(rbinom(500 * 800, 2, 0.3), 500); G = scale(G); ks$grm = tcrossprod(G) / 800   # K 1 = 0
ks$odd = tcrossprod(matrix(rnorm(301 * 301), 301)) / 301                          # odd N
ks$int = { m = outer(id, id, "==") ; storage.mode(m) = "integer"; m }
d = readRDS("/Users/rvialle/suez/bench/data/M_grm.rds"); ks$M_grm = d$sample_kernel
r = lapply(ks, suez:::kernel_contrast_eigen)
saveRDS(r, sprintf("out/eig_%s.rds", a[2]))
