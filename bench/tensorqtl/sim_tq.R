# Simulate the comparison datasets to bench/data/<name>.rds (run from the repository root).
#   Rscript bench/tensorqtl/sim_tq.R [name ...]
source("bench/sim.R")              # make_dataset(); its main block runs only as a script
source("bench/tensorqtl/datasets_tq.R")
dir.create("bench/data", showWarnings = FALSE)
which = commandArgs(TRUE); if (!length(which)) which = names(TQ_DATASETS)
for (nm in which) {
  f = file.path("bench/data", paste0(nm, ".rds"))
  if (file.exists(f)) { cat(nm, ": exists\n"); next }
  d = make_dataset(TQ_DATASETS[[nm]])
  saveRDS(d, f, compress = FALSE)
  cat(nm, ": N =", nrow(d$anno), " genes =", nrow(d$geneloc), " SNPs =", nrow(d$snploc), "\n")
}
