# Round 12+ optimization campaign: shared brief for specialist agents

Repo: /Users/rvialle/suez (R package suez 0.5). Branch `perf-round12`. **Baseline for this campaign = git tag
`perf-baseline-r12`** (suez 0.5 after rounds 1-11, already 2-7x faster than 0.4), installed in
`bench/lib/base12` — never modify or reinstall that library. Older history: `bench/EXPERIMENT_LOG.md`,
`bench/REPORT.md` (read sections 4 and 7 before proposing anything: do not re-propose rejected ideas unless
you have new evidence, and say so explicitly).

## Code
- `R/map_interaction_qtl.R`: per gene (foreach/doMC), cis genotypes -> impute -> filter -> `.Call(suez_fit_gene,
  U, lambda, Uy, UX0, G, C, boot)`; data.frame per gene, `bind_rows` at the end. Set-up once:
  `quantile_normalize`, `kernel_contrast_eigen` (R/lmm.R: Householder + Eigen or LAPACK eigen solver).
- `src/lmm.cpp` (718 lines): the fitting engine. Per gene: eigenvalue-grouped moments, 82-point h grid as a
  GEMM, Schur solves on the shared X0 block, Newton refinement of up to 3 grid peaks per fit, Wald beta/se,
  parametric bootstrap (one y per SNP, then 3 refits). `rotate_all()` = U'G (N x N-1 times N x S) with a
  repeated-row shortcut. Local-ancestry mode: `suez_fit_models` (same machinery, R/local_ancestry.R).
- Step 1 `learn_covariance()` = rstan optimisation (R/learn_covariance.R, src/stan_files) — separate workload.

## Tools
- Datasets in `bench/data/*.rds` (bench/datasets.R, bench/sim.R). Size ladder S (N=200), M (400), L (1000),
  XL (2000), XXL (4000, new); GRM kernels M_grm, L_grm, XL_grm (N=2000, new); M_learn, M_lc; M_num; M_fac3;
  M_eqtl (no interaction); M_cov (10 covariates), M_cov30 (30, new); M_miss (2% missing); M_dense (overlapping
  windows); M_unbal; S_rare; HS; LA (local ancestry); P_many (3000 genes x 10 SNPs, new: per-gene overhead).
- Kernel-level harness: `Sys.setenv(SUEZ_HARNESS_LIB = "bench/lib/base12"); source("bench/proto/harness.R")`
  then `P = prepare("M", n_genes = 30)`, `run_impl(fun, P, reps)`, `compare(a, b)`, `baseline_fit` (= the
  installed base12 .Call). **`bench/proto/r12/p12_base.cpp`** is the base12 engine as a `Rcpp::sourceCpp`
  file exporting `fit_gene_proto(U, lambda, Uy, UX0, G, C, boot)` and `fit_models_proto(...)`: copy it to
  start a C++ prototype so it is drop-in comparable. `bench/proto/r12/p12_timed.cpp` + `profile12.R`: per-phase
  C++ timers (results: `bench/results/r12_profile_base.tsv`).
- End-to-end: `bench/bench.sh <impl> <dataset> <cores> <rep> [step1]` runs the build in `bench/lib/<impl>`
  under /usr/bin/time -l (wall, CPU, max RSS, I/O, process-tree RSS) -> a row in `$RUNS` (default
  bench/results/runs.tsv; use your own file, e.g. RUNS=bench/proto/r12/<role>/runs.tsv). Build a candidate
  package: copy the repo to `bench/proto/r12/<role>/pkg_<name>/` (rsync -a --exclude bench --exclude .git), edit
  there, `R CMD INSTALL --no-docs --library=bench/lib/<role>_<name> bench/proto/r12/<role>/pkg_<name>`.
- R-level profile: `Rscript bench/profile_e2e.R <lib name> <dataset> [n_genes]`.
- End-to-end equivalence: `bench/proto/review/validate.R` (see its header) and `bench/summarize_equivalence.R`.

## Rules
- Work only inside `bench/proto/r12/<your role>/` (and `bench/lib/<role>_*` for builds). Do not edit `src/`,
  `R/`, `tests/` of the repo, do not commit, do not touch other agents' directories or `bench/lib/base12`.
- Correctness vs base12: data log-likelihoods (l0, l_geno, l_interact) |diff| <= 1e-9; beta and se relative
  <= 1e-6 (beta relative to max(|beta|, se)); identical df and NA pattern; bootstrap columns: same R RNG draws
  in the same order, |dl_boot| <= 1e-5 (0.4 itself moves them by 7.6e-7 under a 1e-15 perturbation). Anything
  that changes the statistical model or its outputs beyond that (grid density, number of peaks refined, Newton
  tolerances, skipping fits, a different bootstrap draw order) is a METHODOLOGICAL change: label it, evaluate it
  separately, never present it as a pure speed-up.
- Machine: Apple M5, 10 cores (4 performance + 6 efficiency), 16 GB, macOS, R 4.6 with reference BLAS/LAPACK.
  Other agents benchmark concurrently, so absolute timings are noisy: time your prototype and the baseline
  interleaved in the same process / same script, >= 3 reps, report median ratios. Use small inputs first (S,
  M, M_grm with 20-40 genes); large ones (L, XL, XL_grm, XXL) only for promising prototypes and with few genes.
  Keep at most 2 cores busy yourself unless you are measuring parallel scaling. Final authoritative benchmarks
  are run serially by the coordinator.
- Linux users run suez with OpenBLAS/MKL (see bench/results/linux/): keep that in mind for BLAS-dependent ideas.
- Write findings to `bench/proto/r12/<role>/REPORT.md`: findings, prototypes (file names), correctness
  numbers, speed/RAM ratios with dataset sizes, risks, recommendation. Summarise them in your final message
  (under 400 words) with a ranked list of prototypes.
