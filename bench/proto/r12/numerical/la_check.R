# fit_models (local-ancestry mode) of a prototype vs base12 on captured per-gene inputs of the LA dataset.
#   Rscript bench/proto/r12/numerical/la_check.R <proto.cpp> [genes = 20]
a = commandArgs(TRUE); ng = if (length(a) > 1) as.integer(a[2]) else 20
suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/base12"); library(doMC) }); registerDoMC(1)
CAP = list()
capture_call = function(sym, ...) {
  if (length(CAP) < ng) CAP[[length(CAP) + 1]] <<- list(args = list(...), seed = .Random.seed)
  .Call(sym, ...)
}
f = suez:::la_gene; txt = sub(".Call(suez_fit_models,", "capture_call(suez_fit_models,", deparse(f), fixed = TRUE)
g = eval(parse(text = txt)); environment(g) = asNamespace("suez"); assignInNamespace("la_gene", g, "suez")
d = readRDS("bench/data/LA.rds"); set.seed(42)
d$geneloc = d$geneloc[seq_len(ng), ]
invisible(map_interaction_qtl(d$input, d$genotype, d$geneloc, d$snploc, d$anno, sample_kernel = d$sample_kernel,
                              normalization_approach = "qq", permutation_approach = "boot", cisdist = 1e5, local_ancestry = d$local_ancestry))
env = new.env(); Rcpp::sourceCpp(a[1], env = env)
sym = getNativeSymbolInfo("suez_fit_models", "suez")
run = function(fun) { t0 = proc.time()[["elapsed"]]; r = lapply(CAP, function(c) { assign(".Random.seed", c$seed, envir = globalenv()); do.call(fun, c$args) }); list(r = r, t = proc.time()[["elapsed"]] - t0) }
tb = tp = numeric(3)
for (i in 1:3) { b = run(function(...) .Call(sym, ...)); p = run(env$fit_models_proto); tb[i] = b$t; tp[i] = p$t }
dl = max(sapply(seq_along(CAP), function(i) max(abs(b$r[[i]]$ll - p$r[[i]]$ll), na.rm = TRUE)))
dlb = max(sapply(seq_along(CAP), function(i) suppressWarnings(max(abs(b$r[[i]]$ll_boot - p$r[[i]]$ll_boot), na.rm = TRUE))))
nam = sum(sapply(seq_along(CAP), function(i) sum(is.na(b$r[[i]]$ll) != is.na(p$r[[i]]$ll)) + sum(b$r[[i]]$npar != p$r[[i]]$npar)))
db = max(sapply(seq_along(CAP), function(i) max(mapply(function(x, y, s) max(c(0, abs(x - y) / pmax(abs(x), s)), na.rm = TRUE), b$r[[i]]$beta, p$r[[i]]$beta, b$r[[i]]$se))))
cat(sprintf("LA fit_models: %d genes, base %.3f s, proto %.3f s, ratio %.2f; max |dll| %.1e, |dlb| %.1e, dbeta/max(|b|,se) %.1e, NA/npar mismatches %d\n",
            length(CAP), median(tb), median(tp), median(tb / tp), dl, dlb, db, nam))
