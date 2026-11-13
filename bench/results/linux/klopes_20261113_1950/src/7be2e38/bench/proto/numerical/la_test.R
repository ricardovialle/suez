# fit_models_proto vs the installed suez_fit_models on captured LA inputs: max |diff| per output, NA mismatches, timing
setwd("/Users/rvialle/suez"); suppressPackageStartupMessages(library(suez, lib.loc = "bench/lib/baseline"))
e = new.env(); Rcpp::sourceCpp("bench/proto/numerical/p9_combined.cpp", env = e, cacheDir = "bench/proto/numerical/.cache")
X = readRDS("bench/proto/numerical/la_inputs.rds")
run = function(f) { t0 = proc.time()[[3]]; out = lapply(seq_along(X), function(i) { z = X[[i]]; set.seed(100 + i)
  f(z$U, z$lambda, z$Uy, z$UX0, z$cols, z$models, if (z$boot) 3L else -1L, if (z$boot) 4L else -1L) }); list(out = out, t = proc.time()[[3]] - t0) }
base = function(...) .Call("suez_fit_models", ..., PACKAGE = "suez")
tb = tp = c(); for (r in 1:3) { a = run(base); b = run(e$fit_models_proto); tb[r] = a$t; tp[r] = b$t }
dm = function(get) { x = unlist(lapply(a$out, get)); y = unlist(lapply(b$out, get)); c(max = max(abs(x - y), na.rm = TRUE), na_mis = sum(is.na(x) != is.na(y))) }
print(rbind(ll = dm(function(o) o$ll), npar = dm(function(o) o$npar), ll_boot = dm(function(o) o$ll_boot),
            beta = dm(function(o) unlist(o$beta)), se = dm(function(o) unlist(o$se))))
rb = unlist(lapply(seq_along(a$out), function(i) unlist(a$out[[i]]$beta))); sb = unlist(lapply(a$out, function(o) unlist(o$se)))
pb = unlist(lapply(b$out, function(o) unlist(o$beta))); cat("max |dbeta|/se:", max(abs(rb - pb) / sb, na.rm = TRUE), "\n")
cat("time baseline", median(tb), " proto", median(tp), " speed-up", median(tb / tp), "\n")
