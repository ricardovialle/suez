# Is a sourceCpp prototype timed fairly against the installed base12 .Call? The package is built with
# -DEIGEN_NO_DEBUG (src/Makevars); Rcpp::sourceCpp is not, so Eigen's bounds asserts stay on in p12_base.cpp.
#   Rscript bench/proto/r12/review/flags_check.R [datasets...]
Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12")
source("bench/proto/harness.R"); source("bench/proto/review/interleave.R")
td = tempfile(); dir.create(td)
build = function(name, flags) {
  f = file.path(td, paste0(name, ".cpp")); file.copy("bench/proto/r12/p12_base.cpp", f)
  Sys.setenv(PKG_CPPFLAGS = flags); e = new.env(); Rcpp::sourceCpp(f, env = e, rebuild = TRUE); Sys.unsetenv("PKG_CPPFLAGS")
  e$fit_gene_proto
}
dflt = build("p12_default", "")
nodbg = build("p12_nodebug", "-DEIGEN_NO_DEBUG")
a = commandArgs(TRUE); DS = if (length(a)) a else c("S", "M_grm", "M_cov30")
NG = c(S = 40, M_grm = 30, M_cov30 = 10, L_grm = 8)
for (ds in DS) {
  P = prepare(ds, NG[[ds]])
  cat("==", ds, "\n")
  time_interleaved(list(installed = baseline_fit, sourceCpp_default = dflt, sourceCpp_EIGEN_NO_DEBUG = nodbg), P, reps = 3)
}
