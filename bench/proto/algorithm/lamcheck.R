source("bench/proto/harness.R")
for (ds in c("S","M_fac3","S_rare","M_grm","M_learn","L")) {
  P = prepare(ds, 2)
  l = P$lambda; mx = max(l)
  cat(ds, "n=", length(l), " exact0=", sum(l==0), " tiny(0,1e-8mx]=", sum(l>0 & l<=1e-8*mx),
      " distinct(rel1e-12, ignoring tiny)=", length(unique(signif(l[l>1e-8*mx], 11))), " range tiny:", format(range(l[l>0 & l<=1e-8*mx])), "\n")
  cat("  p0=", ncol(P$UX0), " q=", ncol(P$C), " S/gene=", sapply(P$genes, function(g) ncol(g$G)), "\n")
}
