# suez vs tensorQTL

Scripts behind [doc/suez_vs_tensorqtl.html](https://ricardovialle.github.io/suez/doc/suez_vs_tensorqtl.html): suez and
tensorQTL (CPU and GPU) on the same simulated datasets, run with [pixi](https://pixi.sh) on an otherwise idle machine.

**New round with the current suez** (one command; pulls, builds, tests, re-times, rebuilds the report):

```bash
git pull && mkdir -p bench/lib && bash bench/tensorqtl/rerun.sh 2>&1 | tee bench/lib/rerun.log       # suez only (~30 min)
git pull && mkdir -p bench/lib && bash bench/tensorqtl/rerun.sh all 2>&1 | tee bench/lib/rerun.log   # all, tensorQTL too (~1 h)
```

It stops if the build or `tests/` fail, simulates and exports any missing dataset, keeps the previous outputs in
`bench/tensorqtl/out_<suez commit>_<date>/` (git-ignored), and ends with the speed table and the command to commit
`results/` and `doc/suez_vs_tensorqtl.html`.

`REAL` is the human-cohort shape: 500 individuals x 2 conditions, GRM kernel, 1 gene per 200 kb and a +-1 Mb cis
window (~5,000 SNPs per gene, each SNP in ~10 genes' windows; 150 genes, 734k tests). It runs once, after the
interleaved repetitions, with the rotate-once memory cap at its default (2 GB, `suez_1core` / `suez_8cores`) and at
50 MB (`suez_*_mem5e7`, through `SUEZ_ROTATE_MEM`); tensorQTL uses the same window (`data/REAL/window.txt`).

**Step by step** (the core of `rerun.sh`, from the repository root):

```bash
P() { pixi run --manifest-path bench/tensorqtl/pixi.toml "$@"; }   # function, not a string: works in bash and zsh
P -e r R CMD INSTALL --no-docs --library=bench/lib/cur .        # suez under test
P -e r Rscript bench/tensorqtl/sim_tq.R                          # datasets -> bench/data/
P -e r Rscript bench/tensorqtl/export.R S M L XL XXL G M_grm L_grm XL_grm M_fac3 M_num M_cov M_eqtl M_unbal HS NULLR REAL
bash bench/tensorqtl/run_all.sh                                   # ~1 h; resumes, skipping finished runs
P -e py python bench/tensorqtl/analyze.py                        # results/ and doc/suez_vs_tensorqtl.html
```

| File | Role |
|---|---|
| `pixi.toml` | two environments: `r` (R 4.4, rstan, compilers) and `py` (Python 3.11, CUDA PyTorch, tensorQTL 1.0.10) |
| `datasets_tq.R`, `sim_tq.R` | dataset definitions (most from `bench/datasets.R`) and simulation |
| `export.R` | writes each dataset as tensorQTL inputs (quantile-normalized expression, int8 genotypes) |
| `run_suez.R`, `run_tensorqtl.py` | one timed run each, per-pair p-values and a JSON of timings |
| `run_all.sh` | every run: speed (3 interleaved repetitions), `map_cis` permutations, statistics |
| `rerun.sh` | a new round: build and test suez, missing datasets, back up old outputs, `run_all.sh`, `analyze.py` |
| `analyze.py`, `report_template.html` | summary tables in `results/` and the HTML report |

Timings use monotonic clocks (`/proc/uptime` in R, `perf_counter` in Python): the wall clock jumps under WSL2.
