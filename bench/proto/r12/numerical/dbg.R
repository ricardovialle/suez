Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12")
src = readLines("bench/proto/review/cand/extreme_cases.R"); eval(parse(text = src[grep("^source", src):(grep("^options", src) - 1)]))
P = cases$default_x1e3_opt_x4; g = P$genes[[2]]
for (f in c("bench/proto/r12/numerical/dbg_base.cpp", "bench/proto/r12/numerical/dbg_vec.cpp")) {
  e = new.env(); Rcpp::sourceCpp(f, env = e); e$set_dbg(5L); cat("====", basename(f), "\n"); set.seed(g$seed)
  r = e$fit_gene_proto(P$U, P$lambda, g$Uy, P$UX0, g$G, P$C, FALSE); cat(sprintf("l_geno %.15g\n", r$res[6, 1]))
}
cat(sprintf("exact %.15g\n", exact(cbind(P$UX0, crossprod(P$U, std(g$G[, 6]))), g$Uy, P$lambda)))
