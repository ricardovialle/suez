# quantile_normalize: baseline (apply over rows + t) vs a row loop into a preallocated matrix.
# Rscript bench/proto/memory/qn_test.R <genes> <samples> <which>
a = commandArgs(TRUE); ng = as.integer(a[1]); ns = as.integer(a[2]); which = a[3]
qqnorm_no_plot = function(x) { n = length(x); a = if (n <= 10) 3/8 else 0.5; qnorm((rank(x) - a) / (n + 1 - 2 * a)) }
qn_base = function(input) t(apply(input, 1, qqnorm_no_plot))
qn_loop = function(input) {
  out = matrix(0, nrow(input), ncol(input), dimnames = dimnames(input))
  for (i in seq_len(nrow(input))) out[i, ] = qqnorm_no_plot(input[i, ])
  out
}
rss = function() as.numeric(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE)) / 1024
set.seed(1); x = matrix(rnorm(ng * ns), ng, ns, dimnames = list(paste0("g", 1:ng), paste0("s", 1:ns)))
invisible(gc(reset = TRUE)); r0 = rss()
t = system.time(y <- if (which == "base") qn_base(x) else qn_loop(x))[["elapsed"]]
g = gc()
cat(sprintf("%s %dx%d input %.0f MB: time %.2f s, R heap max %.0f MB (input+output = %.0f MB), RSS %.0f -> %.0f MB, digest %s\n",
  which, ng, ns, object.size(x) / 2^20, t, g[2, ncol(g)], 2 * object.size(x) / 2^20, r0, rss(),
  if (which == "check") identical(qn_base(x), y) else sum(y)))
