suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/base12"); library(doMC) }); registerDoMC(1)
d = readRDS("bench/data/P_many.rds"); g = head(d$geneloc$geneid, 1000)
f0 = suez::map_interaction_qtl
src = deparse(f0); i = grep('cat(gene, nrow(imp_geno)', src, fixed = TRUE); stopifnot(length(i) == 1)
f1 = eval(parse(text = src[-i]))
environment(f1) = asNamespace("suez")
run = function(f) { set.seed(42); t = system.time(r <- f(d$input[g, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno))[["elapsed"]]; list(t = t, r = r) }
for (k in 1:3) { a = run(f0); b = run(f1); cat("base", a$t, "nocat", b$t, "identical", identical(a$r, b$r), "\n") }
