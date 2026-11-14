rss = function() as.numeric(system(sprintf("ps -o rss= -p %d", Sys.getpid()), intern = TRUE)) / 1024
a = commandArgs(TRUE); r0 = rss()
for (p in a) suppressPackageStartupMessages(loadNamespace(p, lib.loc = c(Sys.getenv("LIBX", "bench/lib/base12"), .libPaths())))
invisible(gc()); r1 = rss()
cat(sprintf("%-40s %6.1f MB -> %6.1f MB (+%.1f)\n", paste(a, collapse = "+"), r0, r1, r1 - r0))
