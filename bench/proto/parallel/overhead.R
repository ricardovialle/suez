# Fork / task-dispatch / result-transfer overhead of doMC in a process that has suez and a dataset loaded
# (fork cost grows with the parent's dirty memory). Rscript bench/proto/parallel/overhead.R
suppressPackageStartupMessages({ library(suez, lib.loc = "bench/lib/baseline"); library(doMC); library(foreach) })
d = readRDS("bench/data/M_dense.rds")
tm = function(expr) { t0 = proc.time(); force(expr); t = proc.time() - t0; c(wall = t[["elapsed"]], sys = t[["sys.self"]] + t[["sys.child"]]) }
n = 400
for (cores in c(4, 8)) {
  registerDoMC(cores)
  for (rep in 1:3) {
    s = tm(foreach(i = 1:n) %dopar% i)
    dy = tm(foreach(i = 1:n, .options.multicore = list(preschedule = FALSE)) %dopar% i)
    ch = tm(foreach(i = split(1:n, cut(1:n, 4 * cores, labels = FALSE)), .options.multicore = list(preschedule = FALSE)) %dopar% lapply(i, identity))
    cat(sprintf("cores %d rep %d  n=%d empty tasks: static %.3f s | dynamic (fork/task) %.3f s = %.2f ms/task, sys %.2f s | chunk(4/worker) %.3f s\n",
                cores, rep, n, s[1], dy[1], 1000 * dy[1] / n * cores, dy[2], ch[1]))
  }
}
# result transfer: one M_dense gene result (~1000 SNP rows x 13 columns) serialized back to the parent
r = data.frame(gene = "g", cis_snp = paste0("rs", 1:1000), matrix(rnorm(11000), 1000), stringsAsFactors = FALSE)
cat("serialized size of a 1000-SNP gene result:", length(serialize(r, NULL)) / 1024, "KiB; serialize+unserialize",
    system.time(for (i in 1:100) unserialize(serialize(r, NULL)))[["elapsed"]] * 10, "ms\n")
