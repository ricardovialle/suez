# Kernel-level fits of one installed build (own set-up: the zz sparse basis for the default kernel when the build has it).
#   Rscript integ_run.R <lib> <dataset> <genes> <out.rds>
a = commandArgs(TRUE); lib = a[1]; ds = a[2]; ng = as.integer(a[3]); out = a[4]
Sys.setenv(SUEZ_HARNESS_LIB = lib); source("bench/proto/harness.R")
P = prepare(ds, ng)
d = readRDS(file.path("bench/data", paste0(ds, ".rds")))
if (is.null(d$sample_kernel) && exists("zz_contrast_basis", envir = asNamespace("suez"))) {   # as c1's map_interaction_qtl
  e = suez:::zz_contrast_basis(d$anno$individual); X0 = qr.Q(qr(P$U)) # placeholder, replaced below
  ann = d$anno; one = length(unique(ann$condition)) == 1
  X0 = if (one) matrix(0, nrow(ann), 0) else model.matrix(~ condition, data = ann)
  if (!is.null(d$covariates)) X0 = cbind(X0, d$covariates)
  X0 = suez:::standardize_cols(X0[, apply(X0, 2, sd) > 0, drop = FALSE])
  input = suez:::quantile_normalize(d$input)
  P$U = e$U; P$lambda = e$values; P$UX0 = suez:::crossprod_U(e$U, X0)
  P$genes = lapply(P$genes, function(g) { g$Uy = drop(suez:::crossprod_U(e$U, as.numeric(input[g$gene, ]))); g })
}
r = run_impl(baseline_fit, P); saveRDS(r, out); cat(sprintf("%s %s time %.3f\n", lib, ds, r$time))
