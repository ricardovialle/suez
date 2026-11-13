# Parallel specialist — report (saved by the coordinator from the agent's hand-back)

Recommendation: keep static doMC scheduling; C++ threads only as an opt-in memory saver. All 171 prototype
outputs identical() to baseline 1-core output. Machine shared (load 3-21): ~+-15% noise.

Baseline scaling (median 3 reps): S 1/2/4/8 cores 3.60/2.10/1.36/0.90 s (speed-up 4.0 at 8, eff .50);
M 20.3/11.3/7.22/5.54 s (3.67x); M_dense 1/4/8 54.7/20.4/15.6 s (3.51x). CPU util ~= cores, so loss is extra
CPU time, not idle: efficiency cores are ~2.4x slower per gene (taskpolicy test) -> machine ~6.5 P-core
equivalents; ideal ~5.6x at 8 workers. Gene-size CV 0.02-0.08 in synthetic data.

Scheduling (paired wall ratio baseline/variant; >1 = faster): dynamic (fork per gene) 0.79-1.00 (regression:
14-33 ms per fork), chunk4/chunk8 0.84-1.12 (no gain within noise).

C++ std::thread prototype (grid points and SNP refine loops; draws serial first; bitwise identical): in-process
speed-up 3.1x at 8 threads; end-to-end with 8 CPUs loses to 8 forks (M 9.8 s vs 6.1 s) because of serial
per-gene parts. Memory: 8 forked workers sum to 2.1 GB on M (~220 MB private each: rstan + GC writes break
copy-on-write) vs 0.35 GB for 1 process x 8 threads.

Files: pkg/, prototype.diff, run_variant.sh, scale_baseline.sh, variants.sh, summarize.R, threads_kernel.R,
overhead.R, memsample.sh, pe_speed.R, logs/.
