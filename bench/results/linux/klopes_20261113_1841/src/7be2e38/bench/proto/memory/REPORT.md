# Memory specialist — report (saved by the coordinator from the agent's hand-back)

Breakdown: library(suez) 276 MB / 1.0 s, of which rstan namespace +150-185 MB (QuickJSR stanc heap 76 MB,
ggplot2/loo/...); step 2 never uses rstan. N x N set-up (kernel_contrast_eigen ~10 N^2 doubles; isSymmetric
+138 MB at N = 2000; qr.Q complete +122 MB); U copied per gene by Rcpp::as (32 MB at N = 2000);
quantile_normalize apply+t = input + 6x; result table 125-156 B/row, ~400 MB transient per 1M rows at the end;
checkpoint resume 40 ms/file and not identical (tibble, df_interact double, ~1e-14 rounding).
doMC: GC writes break copy-on-write, ~217 MB private per worker (144 MB with lazy rstan).
Step 1 is the largest RAM user: Stan autodiff (2P + 3 P_adj + 5) N^2 nodes per gradient: S 383 MB, M 663 MB,
L 2.76 GB, XL ~10 GB (extrapolated).

Prototypes (pkg, pkg_lean, pkg_stan; libs bench/lib/mem_*):
- mem_lazy: NAMESPACE without import(rstan); rstan::optimizing in learn_covariance.
- mem_lean: + row-loop quantile_normalize, default kernel as Z Z', kernel freed after eigen, U as Eigen::Map.
- mem_lean2: + column-wise is_symmetric (same all.equal rule), Q via qr.qy (bitwise = qr.Q(complete)[, -1]).
- mem_stan: Stan model with one xflat * s product and one tcrossprod for factors.
mem_lean2 output identical() to baseline on S, M, XL, M_grm, M_learn, M_fac3, M_eqtl, M_num, M_miss, M_cov, LA;
identical at 2/4/8 cores; R CMD check clean (environment warning only).

Measurements (median 3; RSS/footprint MB, 1 core): S 326/308 -> 195/180; M 359/362 -> 230/196;
L 513/516 -> 386/383; XL 704/965 -> 711/721. Load 276 MB/1.0 s -> 129 MB/0.4 s.
Multi-core M total footprint: 1/2/4/8 cores 356/815/1203/2089 -> 223/551/893/1373 (lazy) MB.
Step 1 mem_stan: S 372 -> 326 MB, 2.55 -> 2.02 s; M 686 -> 431 MB, 17.9 -> 9.5 s; L peak 2.77 -> 1.12 GB.
Same log density at the baseline optimum (-10306.2846967856), but the optimizer path differs: cov_suez
1.4e-4 (S) / 3.1e-4 (M) relative; objective equal (S) or marginally better (M).

Risks: future unprefixed rstan calls; regenerating NAMESPACE from an rstantools template; is_symmetric at
exactly the tolerance; mem_stan changes step-1 output ~1e-4.
Recommendation: lazy rstan; mem_lean2; mem_stan once the ~1e-4 step-1 change is accepted; batch genome-wide runs.
