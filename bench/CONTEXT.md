# Optimization project context (shared brief for specialist agents)

Repo: /Users/rvialle/suez (R package, branch `perf-optimization`; baseline = git tag `perf-baseline`,
installed in `bench/lib/baseline` — never modify or reinstall that library).

Hot path: `map_interaction_qtl()` (R/map_interaction_qtl.R) calls `.Call(suez_fit_gene, U, lambda, Uy,
UX0, G, C, boot)` once per gene (src/lmm.cpp). Local-ancestry mode uses `suez_fit_models` (same machinery).
Step 1 `learn_covariance()` is a Stan optimisation (R/learn_covariance.R, src/stan_files).

## Baseline profile (bench/proto/p0_timed.cpp, share of C++ time, 1 core)
| dataset (n contrasts) | rotate U'G | grid (data) | refine (data) | boot sim | boot grids | boot refine |
| S (199)   | 2.4% | 24.8% | 13.4% | 1.6% | 41.7% | 16.2% |
| M_grm (399) | 4.9% | 21.2% | 14.5% | 1.7% | 40.1% | 17.5% |
| L (999)   | 11.2% | 19.1% | 12.6% | 1.6% | 40.1% | 15.4% |
=> grid_values() ~60-65%, refine()/newton ~30%, rotation grows with n. R-side overhead is small for
these sizes. End-to-end 1 core: S 3.3 s, M 18.2 s, L 49.3 s. Loading suez costs ~220 MB RSS (rstan,
StanHeaders) on top of bare R (70 MB): peak RSS on S is 323 MB.

## Tools
- `bench/datasets.R`, `bench/sim.R`: datasets (already generated in bench/data/*.rds; S, M, L, XL size
  ladder; M_grm/M_learn kernels; M_num numeric condition; M_fac3 3 levels; M_eqtl one condition; M_cov 10
  covariates; M_miss missing genotypes; M_dense overlapping cis windows; S_rare rare SNPs; LA local ancestry).
- `bench/proto/harness.R`: `prepare(dataset, n_genes)`, `run_impl(fun, P, reps)`, `compare(a, b)`,
  `baseline_fit` (= installed baseline .Call). `bench/proto/p0_baseline.cpp` is the baseline as a
  `Rcpp::sourceCpp` file exporting `fit_gene_proto(U, lambda, Uy, UX0, G, C, boot)`: copy it to start a
  prototype so that prototypes are drop-in comparable.
- `bench/bench.sh <impl> <dataset> <cores> <rep>`: end-to-end run of the suez build in bench/lib/<impl>
  under /usr/bin/time -l, appends to bench/results/runs.tsv. Install a candidate package build with
  `R CMD INSTALL --no-docs --library=bench/lib/<name> <copy of the package>` (copy the repo to your own
  directory under bench/proto/<role>/pkg first; never edit src/ or R/ of the repo itself).

## Rules
- Work only inside bench/proto/<your role>/. Do not edit src/, R/, tests/, or commit anything.
- Correctness tolerances vs baseline (compare()): log-likelihoods (l0, l_geno, l_interact, bootstrap
  columns) |diff| <= 1e-6; beta/se relative diff <= 1e-6; identical df and NA pattern; bootstrap must
  consume R's RNG in the same order (same draws). Anything that changes the statistical model (grid
  density, number of peaks refined, tolerances, skipping fits) is a METHODOLOGICAL change: label it as
  such, evaluate separately, never present it as a pure speed-up.
- Other agents run benchmarks concurrently on this 10-core machine (4 performance + 6 efficiency cores),
  so absolute timings are noisy: always time your prototype and p0 interleaved in the same process,
  >= 3 reps, report median ratios. Use small inputs (S, M_grm with 20-40 genes) first; L only for
  promising prototypes. Final authoritative benchmarks are run serially by the coordinator.
- Write your findings to bench/proto/<role>/REPORT.md (findings, prototypes with file names, correctness
  numbers, speed/RAM ratios, risks, recommendation) and summarise them in your final message.
