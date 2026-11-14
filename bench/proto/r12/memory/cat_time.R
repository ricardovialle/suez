suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/base12"); library(doMC) }); registerDoMC(1)
d = readRDS("bench/data/P_many.rds"); g = head(d$geneloc$geneid, 1000)
acc = new.env(); acc$t = 0; acc$n = 0
e = new.env(parent = asNamespace("suez"))
e$cat = function(...) { t0 = proc.time()[[3]]; base::cat(...); acc$t = acc$t + proc.time()[[3]] - t0; acc$n = acc$n + 1 }
e$.Call = function(...) { t0 = proc.time()[[3]]; r = base::.Call(...); acc$c = (if (is.null(acc$c)) 0 else acc$c) + proc.time()[[3]] - t0; r }
f = suez::map_interaction_qtl; environment(f) = e
set.seed(42); t = system.time(r <- f(d$input[g, , drop = FALSE], d$genotype, d$geneloc, d$snploc, d$anno))[[3]]
message("wall ", t, " cat total ", acc$t, " over ", acc$n, " calls; .Call total ", acc$c)
