# Round-12 validation suite: one command per candidate build, machine-readable verdict.
#   Rscript bench/proto/r12/review/validate12.R <cand> [--ref base12] [--ds all|quick|S,M_grm,...]
#           [--parts equiv,stat,repro,edge,stress,tests] [--boot exact|dist] [--jobs 2] [--rerun] [--src <pkg dir>]
# <cand>, <ref>: library names in bench/lib (or paths). Output: bench/proto/r12/review/results/v12_<cand>_vs_<ref>.tsv
# (cand ref part dataset check value threshold status detail) and a printed summary ending "VERDICT: PASS|FAIL".
# Parts
#   equiv  end-to-end map_interaction_qtl, candidate vs reference on every dataset (+ S with permutation_approach
#          "none"/"permute", + null runs with permuted genotypes on S, M_grm): rows, NA/Inf pattern, l <= 1e-9,
#          beta/se relative <= 1e-6 (beta relative to max(|beta|, se)), df identical, other numeric columns <= 1e-12;
#          lead SNP per gene, significant sets p < 1e-3/1e-5/1e-8, eGenes (gene p = Bonferroni within gene;
#          Bonferroni across genes and BH 5%), QTL ranks (discordant SNP pairs per gene). Differences inside a 1e-6
#          relative band around a tie or a threshold are listed as INFO "borderline", not failures. Bootstrap columns:
#          --boot exact (default) requires |dl| <= 1e-5 (same draws); --boot dist (candidate that changes the draws,
#          must be declared) gates on the distribution instead (two-sample KS of the bootstrap LR, quantiles).
#   stat   power / FDR vs the simulated truth (identical discoveries expected for a pure speed-up; WARN otherwise),
#          null calibration (lambda_GC, KS) of the candidate's permuted-genotype runs.
#   repro  identical() across 2 runs and across 1/2/4 cores (S, M_dense, M_grm, LA).
#   edge   bench/proto/review/edge_cases.R inputs, candidate and reference .Call out of process (edge12.R).
#   stress random-spectrum and bimodal likelihoods vs an exact optimiser (stress12.R).
#   edge also runs kscale12.R: kernel scaled by 1, 1e6, 1e-6 (tests/model_checks.R) and boundary optima h = 0 / h = 1
#          at scales 1, 1e3, 1e6, vs an exact optimiser over log10 delta in [-25, 25] + both boundaries.
#   tests  every tests/*.R of the candidate's source (--src, else bench/proto/r12/<role>/pkg_<name> for lib <role>_<name>,
#          else the repo's tests/) against the candidate lib; FAIL on any FAIL line or error.
# Reference outputs are cached in bench/proto/r12/review/out/ (delete them if the reference build changes).
suppressPackageStartupMessages(library(parallel))
source("bench/proto/review/validate.R")   # lead_p(), egenes(), add_boot_p(), power_fdr(), calib(), read_tab()
RV = "bench/proto/r12/review"
BENCH_DS = c("S", "M", "L", "XL", "XXL", "M_grm", "L_grm", "XL_grm", "M_learn", "M_lc", "M_num", "M_fac3", "M_eqtl", "M_cov",
             "M_cov30", "M_miss", "M_dense", "M_unbal", "S_rare", "HS", "LA", "P_many")
REV_DS = c("P_many10", "XL_grm_many", "R_real")
QUICK_DS = c("S", "M_grm", "M_learn", "M_fac3", "M_eqtl", "M_cov", "M_miss", "M_unbal", "S_rare", "HS", "LA", "R_real")
TOL = list(ll = 1e-9, boot = 1e-5, rel = 1e-6, other = 1e-12, band = 1e-6)

opt = function(args, name, default) { i = match(paste0("--", name), args); if (is.na(i)) default else args[i + 1] }
a = commandArgs(TRUE)
if (!length(a) || startsWith(a[1], "--")) stop("usage: see the header of bench/proto/r12/review/validate12.R")
cand = a[1]; ref = opt(a, "ref", "base12")
libdir = function(x) if (grepl("/", x)) x else file.path("bench/lib", x)
dsopt = opt(a, "ds", "all")
DS = if (dsopt == "all") c(BENCH_DS, REV_DS) else if (dsopt == "quick") QUICK_DS else strsplit(dsopt, ",")[[1]]
parts = strsplit(opt(a, "parts", "equiv,stat,repro,edge,stress,tests"), ",")[[1]]
boot_mode = opt(a, "boot", "exact"); jobs = as.integer(opt(a, "jobs", "2")); rerun = "--rerun" %in% a
dir.create(file.path(RV, "out"), FALSE); dir.create(file.path(RV, "results"), FALSE)
tsv = file.path(RV, "results", sprintf("v12_%s_vs_%s.tsv", basename(cand), basename(ref)))

ROWS = list()
add = function(part, dataset, check, value, threshold, status, detail = "")
  ROWS[[length(ROWS) + 1]] <<- data.frame(cand = basename(cand), ref = basename(ref), part = part, dataset = dataset, check = check,
    value = if (is.numeric(value)) signif(value, 4) else as.character(value), threshold = as.character(threshold),
    status = status, detail = detail, stringsAsFactors = FALSE)

# ------------------------------------------------------------------ runs (subprocesses, cached reference)
# task: dataset, perm, cores, null seed, tag -> output file
outf = function(lib, ds, perm = "boot", cores = 1, null = NA, tag = "")
  file.path(RV, "out", sprintf("%s_%s_%s_c%d%s%s.rds", basename(lib), ds, perm, cores, if (is.na(null)) "" else paste0("_null", null), tag))
run = function(lib, ds, perm = "boot", cores = 1, null = NA, tag = "", force = FALSE) {
  f = outf(lib, ds, perm, cores, null, tag)
  if (file.exists(f) && !force) return(f)
  log = system2("Rscript", c(file.path(RV, "run12.R"), libdir(lib), ds, cores, f, perm, "Inf", if (!is.na(null)) null),
                stdout = TRUE, stderr = TRUE)
  if (!file.exists(f)) { message("run failed: ", basename(f), "\n", paste(tail(log, 5), collapse = "\n")); return(NA_character_) }
  f
}
pmap = function(tasks, fun) unlist(mclapply(tasks, fun, mc.cores = jobs, mc.preschedule = FALSE))
load_res = function(f) { r = readRDS(f); attr(r, "wall_s") = NULL; as.data.frame(r, stringsAsFactors = FALSE) }

# ------------------------------------------------------------------ equivalence of two result tables
hard_discordant = function(x, y, band) {   # SNP pairs ordered differently by x (ref) and y (cand) beyond the tie band
  o = order(x); x = x[o]; y = y[o]
  k = findInterval(x - band * (1 + abs(x)), x, left.open = TRUE)   # indices i with x_i < x_j - band
  pm = cummax(y); sum(k > 0 & pm[pmax(k, 1)] > y)
}
gene_p = function(r, pc) { l = lead_p(r, pc); data.frame(gene = l$gene, p = pmin(1, l$p * l$m)) }

compare12 = function(A, B, label, boot_mode) {
  part = "equiv"
  ka = paste(A$gene, A$cis_snp); kb = paste(B$gene, B$cis_snp)
  nrow_diff = sum(!ka %in% kb) + sum(!kb %in% ka) + sum(duplicated(ka)) + sum(duplicated(kb))
  add(part, label, "rows", nrow_diff, 0, if (nrow_diff == 0) "PASS" else "FAIL",
      sprintf("ref %d rows / %d genes, cand %d / %d", nrow(A), length(unique(A$gene)), nrow(B), length(unique(B$gene))))
  B = B[match(ka, kb), , drop = FALSE]; keep = !is.na(B$gene); A = A[keep, , drop = FALSE]; B = B[keep, , drop = FALSE]
  miss = setdiff(names(A), names(B)); extra = setdiff(names(B), names(A))
  add(part, label, "columns", length(miss) + length(extra), 0, if (length(miss) + length(extra)) "FAIL" else "PASS", paste(c(miss, extra), collapse = " "))
  num = intersect(names(A)[sapply(A, is.numeric)], names(B))
  for (cn in num) {
    x = A[[cn]]; y = B[[cn]]
    pat = sum(is.na(x) != is.na(y) | (is.infinite(x) != is.infinite(y)) | (is.infinite(x) & is.infinite(y) & x != y), na.rm = TRUE)
    ok = is.finite(x) & is.finite(y); d = abs(x - y)[ok]
    boot = grepl("^l_boot", cn)
    kind = if (boot) "boot" else if (cn == "l0" || grepl("^l_", cn)) "ll" else if (grepl("^beta", cn)) "beta" else if (grepl("^se_", cn)) "se" else
      if (grepl("^df_", cn)) "df" else if (grepl("^p_", cn)) "p" else "other"
    if (kind %in% c("beta", "se")) {
      den = pmax(abs(x), abs(y))[ok]
      if (kind == "beta") { s = A[[sub("^beta", "se", cn)]]; if (!is.null(s)) den = pmax(den, s[ok], na.rm = TRUE) }
      d = ifelse(d == 0, 0, d / den)
    }
    mx = if (length(d)) max(d) else 0
    thr = switch(kind, ll = TOL$ll, boot = TOL$boot, beta = , se = TOL$rel, df = 0, p = NA, other = TOL$other)
    st = if (pat > 0) "FAIL" else if (is.na(thr)) "INFO" else if (mx <= thr) "PASS" else if (boot && boot_mode == "dist") "INFO" else "FAIL"
    add(part, label, paste0(cn, if (kind == "p") " max|dlog10p|" else if (kind %in% c("beta", "se")) " max rel" else " max abs"),
        if (kind == "p") { lx = -log10(x[ok]); ly = -log10(y[ok]); f = is.finite(lx) & is.finite(ly); if (any(f)) max(abs(lx - ly)[f]) else 0 } else mx,
        if (is.na(thr)) "-" else thr, st, sprintf("n=%d na/inf mismatches=%d identical=%d", sum(ok), pat, sum(x[ok] == y[ok])))
  }
  # decisions derived from p-values; only data p-values are gated when the bootstrap draws may differ
  A = add_boot_p(A); B = add_boot_p(B)
  pcs = intersect(c("p_geno", "p_interact", "p_het", "p_interact_het", "p_interact_boot"), intersect(names(A), names(B)))
  for (pc in pcs) {
    gate = !(pc == "p_interact_boot" && boot_mode == "dist")
    if (all(is.na(A[[pc]]))) next
    fail = if (gate) "FAIL" else "INFO"
    x = A[[pc]]; y = B[[pc]]
    # lead SNP per gene; a different lead whose reference p is within the band of the reference lead is a tie
    la = lead_p(A, pc); lb = lead_p(B, pc); lb = lb[match(la$gene, lb$gene), ]
    diff = which(la$snp != lb$snp | is.na(lb$snp))
    pa = setNames(x, ka[keep])
    tie = vapply(diff, function(i) { pr = pa[paste(la$gene[i], lb$snp[i])]
      !is.na(pr) && abs(log10(pr) - log10(la$p[i])) <= TOL$band * (1 + abs(log10(la$p[i]))) }, TRUE)
    nh = sum(!tie)
    add(part, label, paste(pc, "lead SNP per gene"), nh, 0, if (nh) fail else "PASS",
        sprintf("%d genes, %d different leads (%d ties within the band)", nrow(la), length(diff), sum(tie)))
    # significant sets
    for (t in c(1e-3, 1e-5, 1e-8)) {
      fl = which(xor(x < t, y < t)); bd = abs(log10(x[fl]) - log10(t)) <= TOL$band * (1 + abs(log10(t)))
      add(part, label, sprintf("%s sig p<%g", pc, t), sum(!bd), 0, if (sum(!bd)) fail else "PASS",
          sprintf("ref %d cand %d flips %d (borderline %d)", sum(x < t, na.rm = TRUE), sum(y < t, na.rm = TRUE), length(fl), sum(bd)))
    }
    # eGenes: Bonferroni within gene, then Bonferroni across genes or BH 5%
    ga = gene_p(A, pc); gb = gene_p(B, pc); gb = gb[match(ga$gene, gb$gene), ]
    for (rule in c("Bonferroni", "BH")) {
      qa = if (rule == "BH") p.adjust(ga$p, "BH") else pmin(1, ga$p * nrow(ga)); qb = if (rule == "BH") p.adjust(gb$p, "BH") else pmin(1, gb$p * nrow(gb))
      fl = which(xor(qa < 0.05, qb < 0.05)); bd = abs(log10(qa[fl]) - log10(0.05)) <= TOL$band * 2
      add(part, label, sprintf("%s eGenes %s 5%%", pc, rule), sum(!bd), 0, if (sum(!bd)) fail else "PASS",
          sprintf("ref %d cand %d flips %d (borderline %d)", sum(qa < 0.05, na.rm = TRUE), sum(qb < 0.05, na.rm = TRUE), length(fl), sum(bd)))
    }
    # rank agreement within genes
    idx = split(seq_along(x), A$gene)
    hd = sum(vapply(idx, function(i) { ok = !is.na(x[i]) & !is.na(y[i]); if (sum(ok) < 2) 0L else
      as.integer(hard_discordant(-log10(pmax(x[i][ok], 1e-300)), -log10(pmax(y[i][ok], 1e-300)), TOL$band)) }, 0L))
    sp = suppressWarnings(vapply(idx, function(i) { ok = !is.na(x[i]) & !is.na(y[i])
      if (sum(ok) < 3) NA_real_ else cor(x[i][ok], y[i][ok], method = "spearman") }, 0))
    add(part, label, paste(pc, "rank: discordant SNP pairs"), hd, 0, if (hd) fail else "PASS",
        sprintf("Spearman per gene min %.6f median %.6f", suppressWarnings(min(sp, na.rm = TRUE)), median(sp, na.rm = TRUE)))
  }
  # bootstrap distribution (always reported; gates only with --boot dist)
  lr = function(r) { if (all(c("l_boot_interact", "l_boot_geno") %in% names(r))) 2 * (r$l_boot_interact - r$l_boot_geno) else
    if (all(c("l_boot_interact", "l_boot_interact_null") %in% names(r))) 2 * (r$l_boot_interact - r$l_boot_interact_null) else NULL }
  for (what in c("boot LR", "l_boot_geno - l_geno")) {
    u = if (what == "boot LR") lr(A) else if (!is.null(A$l_boot_geno)) A$l_boot_geno - A$l_geno
    v = if (what == "boot LR") lr(B) else if (!is.null(B$l_boot_geno)) B$l_boot_geno - B$l_geno
    if (is.null(u) || sum(is.finite(u)) < 20) next
    u = u[is.finite(u)]; v = v[is.finite(v)]
    ks = suppressWarnings(ks.test(u, v)$p.value); q = c(0.5, 0.9, 0.99, 0.999)
    st = if (boot_mode != "dist") "INFO" else if (ks < 1e-3) "FAIL" else if (ks < 0.01) "WARN" else "PASS"
    add(part, label, paste(what, "distribution KS p"), ks, if (boot_mode == "dist") ">= 0.01" else "-", st,
        paste(sprintf("q%g %.4g/%.4g", q, quantile(u, q), quantile(v, q)), collapse = " "))
  }
}

# ------------------------------------------------------------------ parts
t_start = proc.time()[["elapsed"]]
if ("equiv" %in% parts || "stat" %in% parts) {
  tasks = c(lapply(DS, function(d) list(ds = d)),
            if ("S" %in% DS) list(list(ds = "S", perm = "none"), list(ds = "S", perm = "permute")),
            lapply(intersect(c("S", "M_grm"), DS), function(d) list(ds = d, null = 1)))
  lab = function(t) paste0(t$ds, if (!is.null(t$perm)) paste0("@", t$perm), if (!is.null(t$null)) "#null")
  go = function(lib, force, tag) pmap(tasks, function(t) run(lib, t$ds, if (is.null(t$perm)) "boot" else t$perm, 1,
                                                             if (is.null(t$null)) NA else t$null, tag = tag, force = force))
  fr = go(ref, rerun, ""); fc = go(cand, TRUE, "_cand")   # candidate always rerun (its build may have changed)
  for (i in seq_along(tasks)) {
    L = lab(tasks[[i]])
    if (is.na(fr[i]) || is.na(fc[i])) { add("equiv", L, "run", if (is.na(fc[i])) "cand failed" else "ref failed", "-", "FAIL"); next }
    A = load_res(fr[i]); B = load_res(fc[i])
    wr = attr(readRDS(fr[i]), "wall_s"); wc = attr(readRDS(fc[i]), "wall_s")
    add("equiv", L, "wall_s cand/ref (indicative, not a benchmark)", wc / wr, "-", "INFO", sprintf("%.2f / %.2f s", wc, wr))
    if ("equiv" %in% parts) compare12(A, B, L, boot_mode)
    if ("stat" %in% parts) {
      ds = tasks[[i]]$ds; f = file.path("bench/data", paste0(ds, ".rds")); if (!file.exists(f)) f = file.path(RV, "data", paste0(ds, ".rds"))
      if (is.null(tasks[[i]]$null) && is.null(tasks[[i]]$perm)) {
        truth = readRDS(f)$truth
        if (!is.null(truth)) {
          pa = power_fdr(A, truth); pb = power_fdr(B, truth)
          if (boot_mode == "dist") { pa = pa[pa$test != "p_interact_boot", ]; pb = pb[pb$test != "p_interact_boot", ] }
          same = identical(pa$discoveries, pb$discoveries) && isTRUE(all.equal(pa$power, pb$power)) && isTRUE(all.equal(pa$fdr, pb$fdr))
          add("stat", L, "power/FDR vs truth", sprintf("%d/%d rules identical", sum(pa$discoveries == pb$discoveries & pa$power == pb$power), nrow(pa)), "identical",
              if (same) "PASS" else "WARN", paste(sprintf("%s %s: power %.3f/%.3f fdr %.3f/%.3f", pa$test, pa$rule, pa$power, pb$power, pa$fdr, pb$fdr), collapse = "; "))
        }
      }
      if (!is.null(tasks[[i]]$null)) {
        ca = capture.output(xa <- calib(A)); cb = capture.output(xb <- calib(B))
        for (j in seq_len(nrow(xb))) {   # 1 SNP/gene rows: ~150 genes, lambda's sd ~0.15, so only the cand-ref gap gates them
          lg = xb$lambda_gc[j]; dg = abs(lg - xa$lambda_gc[j])
          st = if (dg > 0.05) "FAIL" else if (dg > 0.02 || (!grepl("1 SNP/gene", xb$test[j]) && (lg < 0.9 || lg > 1.1))) "WARN" else "PASS"
          add("stat", L, paste("null lambda_GC", xb$test[j]), lg, "|cand-ref| <= 0.02 (FAIL > 0.05); all-SNP/boot rows 0.9-1.1", st,
              sprintf("ref %.3f; cand KS p %.3g (ref %.3g)", xa$lambda_gc[j], xb$ks_p[j], xa$ks_p[j]))
        }
      }
    }
  }
}
if ("repro" %in% parts) {
  for (d in c("S", "M_dense", "M_grm", "LA")) {   # always these (small); 4-core runs are short bursts
    f1 = run(cand, d, tag = "_cand", force = !(d %in% DS && "equiv" %in% parts))   # the equiv run of this invocation, or a fresh one
    others = list(c(1, "_cand_rep2"), c(2, "_cand"), c(4, "_cand"))
    fo = vapply(others, function(o) run(cand, d, cores = as.integer(o[1]), tag = o[2], force = TRUE), "")
    r1 = load_res(f1)
    for (k in seq_along(others)) {
      same = !is.na(fo[k]) && identical(r1, load_res(fo[k]))
      add("repro", d, sprintf("identical(1 core, %s)", if (others[[k]][2] == "_cand_rep2") "1 core rerun" else paste(others[[k]][1], "cores")),
          same, TRUE, if (same) "PASS" else "FAIL")
    }
  }
}
ext = function(script, part, extra = character(0)) {
  f = file.path(RV, "results", sprintf("%s_%s.tsv", part, basename(cand)))
  log = system2("Rscript", c(file.path(RV, script), libdir(cand), libdir(ref), f, extra), stdout = TRUE, stderr = TRUE)
  if (!file.exists(f)) { add(part, "-", "run", "script failed", "-", "FAIL", paste(tail(log, 3), collapse = " | ")); return() }
  x = read.delim(f, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(x))) add(part, x$dataset[i], x$check[i], x$value[i], x$threshold[i], x$status[i], x$detail[i])
  invisible(file.remove(f))
}
if ("edge" %in% parts) { ext("edge12.R", "edge"); ext("kscale12.R", "edge") }
if ("tests" %in% parts) {   # every tests/*.R of the candidate's own source against the candidate lib
  guess = sub("^([^_]+)_(.+)$", "bench/proto/r12/\\1/pkg_\\2", basename(cand))
  src = opt(a, "src", if (dir.exists(file.path(guess, "tests"))) guess else ".")
  tf = list.files(file.path(src, "tests"), "\\.R$", full.names = TRUE)
  if (!length(tf)) add("tests", src, "test files", 0, ">= 1", "FAIL", "no tests/*.R")
  for (f in tf) {
    log = suppressWarnings(system2("Rscript", c("-e", shQuote(sprintf("suppressPackageStartupMessages(library(suez, lib.loc = '%s')); source('%s')", libdir(cand), f))),
                                   stdout = TRUE, stderr = TRUE))
    st = attr(log, "status"); fl = grep("^FAIL", log, value = TRUE); np = length(grep("^PASS", log))
    err = !is.null(st) && st != 0
    add("tests", basename(f), "PASS / FAIL lines, error", sprintf("%d / %d%s", np, length(fl), if (err) ", error" else ""), "0 FAIL, no error",
        if (length(fl) || err || np == 0) "FAIL" else "PASS",
        paste(c(sprintf("source %s", src), substr(fl, 1, 150), if (err) substr(tail(log, 2), 1, 150)), collapse = " | "))
  }
}
if ("stress" %in% parts) ext("stress12.R", "stress")

out = do.call(rbind, ROWS)
write.table(out, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s vs %s: %d checks, %s (%.0f s) -> %s\n", basename(cand), basename(ref), nrow(out),
            paste(names(table(out$status)), table(out$status), sep = "=", collapse = " "), proc.time()[["elapsed"]] - t_start, tsv))
bad = out[out$status %in% c("FAIL", "WARN"), ]
if (nrow(bad)) { op = options(width = 250); print(bad[, c("part", "dataset", "check", "value", "threshold", "status", "detail")], row.names = FALSE); options(op) }
cat("VERDICT:", if (any(out$status == "FAIL")) "FAIL" else if (any(out$status == "WARN")) "PASS with warnings" else "PASS", "\n")
